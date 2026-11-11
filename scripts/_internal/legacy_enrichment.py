"""Restore legacy display fields and enrich human protein metadata, with an audit.

Default: offline, JSON/CSV outputs. --online enables the public UniProt REST API;
--xlsx additionally exports a workbook. Original files are never overwritten.
Based on paper_test/download_sequence*.py, download_localization_5.py and
convert_hgvs.py; replaces first-hit selection and incomplete indel conversion.
API reference: https://www.uniprot.org/help/api_queries
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERSION = '1.0.2'
SPEC = importlib.util.spec_from_file_location('legacy_normalizer', ROOT/'scripts/_internal/record_cleaning.py')
N = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(N)
BEFORE, AFTER = N.BASE_COLUMNS[6:8]
ENRICH_FIELDS = ['Accession_new', 'Reviewed', 'Wild Sequence_new', 'Subcellular_Location_Values']
MISSING = {'', 'unknown', 'not_reported', 'not reported', 'null', 'n/a', 'na', '*', 'not found', 'undetermined'}
ACCESSION = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[1-9]\d*)?')


def text(value):
    return N.text(value).strip()


def missing(value):
    return text(value).casefold() in MISSING


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=path.parent, suffix='.tmp')
    try:
        with os.fdopen(handle, 'w', encoding='utf-8') as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def load_rows(path, sheet=None):
    """Load a flat legacy table; never coerce PMID or stable IDs to numbers."""
    if path.suffix.lower() == '.csv':
        with path.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            headers, rows = reader.fieldnames, list(reader)
    elif path.suffix.lower() == '.json':
        payload = read_json(path)
        rows = payload.get('records', payload) if isinstance(payload, dict) else payload
        if not isinstance(rows, list) or not rows or not all(isinstance(r, dict) for r in rows):
            raise ValueError('JSON input must contain an array of legacy-column row objects.')
        headers = list(dict.fromkeys(h for row in rows for h in row))
    else:
        from openpyxl import load_workbook
        # Own the file handle explicitly: partially consumed read-only worksheet
        # iterators can otherwise keep a Windows file locked after book.close().
        with path.open('rb') as source_stream:
            book = load_workbook(source_stream, read_only=True, data_only=True)
            try:
                if sheet:
                    chosen = book[sheet]
                elif '完整记录_26列' in book.sheetnames:
                    chosen = book['完整记录_26列']
                else:
                    candidates = [s for s in book if {'id','PMID','Gene Name','Mutation Information'} <= set(next(s.values, ()))]
                    chosen = candidates[0] if len(candidates)==1 else None
                if chosen is None: raise ValueError('Select the data sheet using --sheet.')
                iterator = chosen.values
                raw_headers = list(next(iterator))
                headers = [h for h in raw_headers if h is not None]
                rows = [{h:r[i] if i<len(r) else None for i,h in enumerate(raw_headers) if h is not None}
                        for r in iterator if any(v is not None for v in r)]
            finally:
                book.close()
    if not headers or len(set(headers)) != len(headers): raise ValueError('Duplicate or empty headers.')
    required = {'id','PMID','Gene Name','Mutation Type','Mutation Information',BEFORE,AFTER,'Mislocalization Phenotype'}
    if not required <= set(headers): raise ValueError('Missing columns: '+', '.join(sorted(required-set(headers))))
    seen = set()
    for row in rows:
        row['id'], row['PMID'] = text(row['id']), text(row['PMID'])
        if not row['id'] or row['id'] in seen or not re.fullmatch(r'[1-9]\d*',row['PMID']):
            raise ValueError('Empty/duplicate ID or invalid PMID: '+row['id'])
        seen.add(row['id'])
    for h in ENRICH_FIELDS:
        if h not in headers: headers.append(h)
    for row in rows:
        for h in headers: row.setdefault(h,None)
    return headers, rows


class HistoricalSources:
    """Bind stable IDs to explicit JSON locators; never match by gene alone."""
    def __init__(self, audit_path=None, original_path=None, dedup_path=None):
        self.records, self.locators, self.hashes = {}, {}, {}
        self.groups = defaultdict(list)
        self.audit_path, self.original_path = audit_path, original_path
        if dedup_path:
            self.hashes[str(dedup_path)] = digest(dedup_path)
            mapping = read_json(dedup_path).get('canonical')
            if not isinstance(mapping, dict): raise ValueError('Dedup map requires a canonical ID mapping.')
            for rid, target in mapping.items(): self.groups[target].append(rid)

    def load(self, rows):
        wanted = {rid for row in rows for rid in self.groups.get(row['id'], [row['id']])}
        if self.audit_path:
            self.hashes[str(self.audit_path)] = digest(self.audit_path)
            with self.audit_path.open(encoding='utf-8-sig', newline='') as stream:
                entries = [r for r in csv.DictReader(stream) if r['id'] in wanted]
            grouped = defaultdict(list)
            for entry in entries:
                path = Path(entry['原始抽取文件'])
                grouped[path if path.is_absolute() else ROOT/path].append(entry)
            for path, group in grouped.items():
                self.hashes[str(path)] = digest(path)
                payload = read_json(path)
                for item in group:
                    locator = item['JSON定位']
                    match = re.fullmatch(r'records\[(\d+)\]|\[(\d+)\]\.original_record',locator)
                    if not match: raise ValueError('Unsupported source locator: '+locator)
                    record = payload['records'][int(match[1])] if match[1] is not None else payload[int(match[2])]['original_record']
                    if text(record.get('pmid')) != text(item['PMID']): raise ValueError('Source PMID mismatch: '+item['id'])
                    self.records[item['id']] = record
                    self.locators[item['id']] = str(path)+'#'+locator
        elif self.original_path:
            self.hashes[str(self.original_path)] = digest(self.original_path)
            payload = read_json(self.original_path)
            counts = Counter()
            for i, record in enumerate(payload['records']):
                pmid = text(record.get('pmid')); counts[pmid] += 1
                rid = f'{pmid}:{counts[pmid]:04d}'
                if rid in wanted:
                    self.records[rid] = record
                    self.locators[rid] = str(self.original_path)+f'#records[{i}]'
        for row in rows:
            for rid in self.groups.get(row['id'], [row['id']]):
                record = self.records.get(rid)
                if record and text(record.get('pmid')) != row['PMID']:
                    raise ValueError('Cross-paper source binding: '+rid)

    def members(self, rid):
        return self.groups.get(rid, [rid])


def location_key(value):
    return tuple(t.strip().casefold() for t in text(value).split(';'))


def restore_other(current, raw_values):
    """Only restore explicit historical Other-X tokens with identical other tokens."""
    if not any(t.strip().casefold()=='other' for t in text(current).split(';')):
        return current, 'NOT_NEEDED'
    proposals = []
    for raw in raw_values:
        if missing(raw): return current, 'SOURCE_DESCRIPTION_MISSING'
        collapsed, restored = [], []
        for token in text(raw).split(';'):
            token = token.strip()
            match = re.fullmatch(r'Other\s*[-–—:]\s*(.+)',token,re.I)
            if match and not missing(match[1]) and match[1].strip().casefold()!='other':
                collapsed.append('Other'); restored.append('Other-'+match[1].strip())
            else:
                collapsed.append(token); restored.append(token)
        if location_key(';'.join(collapsed)) != location_key(current):
            return current, 'RAW_STRUCTURE_DIFFERENT'
        if any(token.casefold()=='other' for token in restored):
            return current, 'SOURCE_DESCRIPTION_INCOMPLETE'
        proposals.append('; '.join(restored))
    if not proposals: return current, 'SOURCE_MISSING'
    if len({location_key(p) for p in proposals}) != 1:
        return current, 'MERGED_SOURCE_DETAILS_DIFFER'
    return proposals[0], 'RESTORED_FROM_RAW'


def normalize_variant(value, mutation_type):
    """Formatting only; no inferred protein consequence, position or cis/trans phase."""
    original = text(value)
    normalized, reason = N.normalize_mutation(original, mutation_type)
    if N.P_RE.fullmatch(text(normalized)) or normalized in {'p.=','p.(=)'}:
        return normalized, 'SINGLE_PROTEIN_FORMAT', reason
    # Complete one-letter indel endpoints and explicitly stated inserted residues.
    token = original.removeprefix('p.')
    predicted = token.startswith('(') and token.endswith(')')
    if predicted: token = token[1:-1]
    match = re.fullmatch(r'([ARNDCEQGHILKMFPSTWYV])([1-9]\d*)(?:_([ARNDCEQGHILKMFPSTWYV])([1-9]\d*))?(delins|ins|del|dup)([ARNDCEQGHILKMFPSTWYV]*)',token)
    if match:
        left,pos,right,end,operation,inserted = match.groups()
        valid = (operation in {'del','dup'} and not inserted) or (operation=='delins' and bool(inserted)) or (operation=='ins' and bool(inserted) and right and int(end)==int(pos)+1)
        valid = valid and (not end or int(end)>int(pos))
        if valid:
            body = N.AA[left]+pos+(('_'+N.AA[right]+end) if right else '')+operation+''.join(N.AA[c] for c in inserted)
            candidate = 'p.'+('('+body+')' if predicted else body)
            if N.P_RE.fullmatch(candidate): return candidate, 'SINGLE_PROTEIN_FORMAT', 'Explicit one-letter indel notation expanded; reference sequence unverified'
    tokens = list(dict.fromkeys(m[0] for m in N.P_RE.finditer(original)))
    category = 'CONTEXT_PRESERVED' if len(tokens)==1 else 'MULTIPLE_VARIANTS_PRESERVED' if len(tokens)>1 else 'OTHER_NOTATION_PRESERVED'
    return original, category, reason


class UniProtClient:
    """Persistent request cache; exact human identities only; errors are not misses."""
    def __init__(self, cache_dir, online=False, max_requests=None, timeout=30, retries=2,
                 interval=0.4, max_age_days=30, include_unreviewed=False, session=None):
        import requests
        self.cache_dir = Path(cache_dir)
        self.online, self.max_requests = online, max_requests
        self.timeout, self.retries, self.interval = timeout, retries, interval
        self.max_age_days, self.include_unreviewed = max_age_days, include_unreviewed
        self.session = session or requests.Session()
        self.calls, self.cache_hits, self.last_request, self.failure_streak = 0, 0, 0, 0
        self.memo = {}

    def request(self, endpoint, params=None):
        import requests
        url = 'https://rest.uniprot.org/uniprotkb/'+endpoint
        identity = {'url':url,'params':params or {}}
        key = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
        if key in self.memo: return self.memo[key]
        path = self.cache_dir/(key+'.json')
        if path.exists():
            cached = read_json(path)
            if cached.get('schema')!='uniprot-request-cache-v1' or cached.get('request')!=identity:
                raise ValueError('Cache identity/schema mismatch: '+str(path))
            age = (datetime.now(timezone.utc)-datetime.fromisoformat(cached['fetched_at'])).total_seconds()/86400
            if age <= self.max_age_days:
                self.cache_hits += 1; self.memo[key] = cached
                return cached
            if not self.online: return {'status':'STALE_CACHE','request':identity}
        if not self.online: return {'status':'OFFLINE_CACHE_MISS','request':identity}
        if self.failure_streak>=3: return {'status':'DEFERRED_AFTER_FAILURES','request':identity}
        for attempt in range(self.retries+1):
            if self.max_requests is not None and self.calls>=self.max_requests:
                return {'status':'REQUEST_BUDGET_EXHAUSTED','request':identity}
            time.sleep(max(0,self.interval-(time.monotonic()-self.last_request)))
            self.last_request = time.monotonic(); self.calls += 1
            try:
                response = self.session.get(url,params=params,timeout=(5,self.timeout),headers={'Accept':'application/json','User-Agent':'MutProMisLoc-enrichment/'+VERSION})
                if response.status_code==429 or response.status_code>=500:
                    if attempt<self.retries:
                        delay = response.headers.get('Retry-After','')
                        time.sleep(min(30,float(delay) if delay.isdigit() else 2**attempt)); continue
                if response.status_code==404:
                    result = {'status':'NOT_FOUND','payload':None}
                else:
                    response.raise_for_status()
                    payload = response.json()
                    if not isinstance(payload,dict): raise ValueError('Expected JSON object')
                    result = {'status':'OK','payload':payload}
                result.update(schema='uniprot-request-cache-v1',request=identity,
                              fetched_at=datetime.now(timezone.utc).isoformat(),
                              release=response.headers.get('X-UniProt-Release'),response_url=response.url,
                              total=response.headers.get('X-Total-Results'),next_page=bool(response.links.get('next')))
                atomic_json(path,result); self.memo[key]=result
                self.failure_streak=0
                return result
            except (requests.RequestException, ValueError) as exc:
                if attempt==self.retries:
                    self.failure_streak+=1
                    result={'status':'REQUEST_FAILED','request':identity,'error':type(exc).__name__+': '+str(exc)}
                    self.memo[key]=result
                    return result
                time.sleep(min(4,2**attempt))

    def resolve(self, gene, accession=None):
        name = N.uniprot_query_name(gene)
        if not name: return {'status':'NAME_REQUIRES_REVIEW','matches':[]}
        if accession and ACCESSION.fullmatch(text(accession)):
            accession = text(accession)
            if '-' in accession: return {'status':'ISOFORM_ACCESSION_REQUIRES_REVIEW','matches':[]}
            response = self.request(accession+'.json')
            if response['status']!='OK': return {**response,'matches':[]}
            entry = response['payload']
            identity_ok = accession==entry.get('primaryAccession') or accession in entry.get('secondaryAccessions',[])
            exact = name.casefold() in {text(v).casefold() for v in N.HumanUniProtLookup.exact_names(entry)}
            if not identity_ok or entry.get('organism',{}).get('taxonId')!=9606 or not exact:
                return {**response,'status':'ACCESSION_IDENTITY_CONFLICT','matches':[]}
            return {**response,'status':'UNIQUE','matches':[entry]}
        for reviewed in ([True,False] if self.include_unreviewed else [True]):
            query = f'organism_id:9606 AND reviewed:{str(reviewed).lower()} AND (gene_exact:"{name}" OR protein_name:"{name}")'
            response = self.request('search',{'query':query,'format':'json','size':100})
            if response['status']!='OK': return {**response,'matches':[]}
            entries = response['payload'].get('results')
            if not isinstance(entries,list): return {**response,'status':'INVALID_RESPONSE','matches':[]}
            total=text(response.get('total'))
            if total and not total.isdigit(): return {**response,'status':'INVALID_RESPONSE','matches':[]}
            if response.get('next_page') or int(total or len(entries))>len(entries) or len(entries)>=100:
                return {**response,'status':'TRUNCATED_SEARCH','matches':[]}
            matches = {e['primaryAccession']:e for e in entries if e.get('primaryAccession') and e.get('organism',{}).get('taxonId')==9606
                       and name.casefold() in {text(v).casefold() for v in N.HumanUniProtLookup.exact_names(e)}}
            if matches: return {**response,'status':'UNIQUE' if len(matches)==1 else 'AMBIGUOUS','matches':list(matches.values())}
        return {**response,'status':'NOT_FOUND_EXACT','matches':[]}

    def close(self): self.session.close()


def sequence_check(variant, sequence):
    """Check stated reference residues only; never renumber or rewrite variants."""
    value=text(variant)
    if not N.P_RE.fullmatch(value):
        prefix=N.P_RE.match(value)
        suffix=value[prefix.end():] if prefix else ''
        if not prefix or not re.fullmatch(r'\s+\([^()]+\)',suffix) or re.search(r'\bp\.',suffix):
            return 'NOT_CHECKED_COMPLEX_NOTATION'
        value=prefix[0]
    positions = re.findall(r'('+'|'.join(N.AA.values())+r')([1-9]\d*)',value)
    inverse = {v:k for k,v in N.AA.items()}
    if not positions: return 'NOT_CHECKED'
    for residue, position in positions:
        p = int(position)
        if p>len(sequence): return 'POSITION_OUT_OF_RANGE'
        if sequence[p-1]!=inverse[residue]: return 'REFERENCE_RESIDUE_MISMATCH'
    return 'STATED_REFERENCE_RESIDUES_MATCH'


def database_location_applicable(row):
    """Only an explicit no-change label; equal compartments are not evidence."""
    phenotype = re.sub(r'[\s_-]+', '', text(row.get('Mislocalization Phenotype')).casefold())
    provenance = row.get('Outcome_Provenance')
    return phenotype == 'nochange' and (missing(provenance) or text(provenance) == 'EXPERIMENTAL_NO_CHANGE')


def entry_values(entry, row=None):
    entry_type = entry.get('entryType','')
    reviewed = 'reviewed' if entry_type.startswith('UniProtKB reviewed') else 'unreviewed' if entry_type.startswith('UniProtKB unreviewed') else None
    locations = []
    for comment in entry.get('comments',[]):
        if comment.get('commentType')!='SUBCELLULAR LOCATION' or comment.get('molecule'): continue
        for location in comment.get('subcellularLocations',[]):
            value = location.get('location',{}).get('value')
            if value and value not in locations: locations.append(value)
    return {'Accession_new':entry.get('primaryAccession'),'Reviewed':reviewed,
            'Wild Sequence_new':entry.get('sequence',{}).get('value'),
            'Subcellular_Location_Values':('; '.join(locations) or None)
            if row is None or database_location_applicable(row) else None}


def process(rows, sources, client):
    output, changes, audit, source_details = [], [], [], []
    for i, original in enumerate(rows,1):
        row = dict(original); rid = row['id']; members = sources.members(rid)
        raw_records = [sources.records.get(member,{}) for member in members]
        notes = []
        def update(field,value,rule):
            if row.get(field)!=value:
                changes.append({'id':rid,'PMID':row['PMID'],'field':field,'before':row.get(field),'after':value,'rule':rule})
                row[field]=value
        if not database_location_applicable(row):
            update('Subcellular_Location_Values', None, 'DATABASE_LOCATION_ONLY_FOR_NO_CHANGE')
        for field,key in [(BEFORE,'normal_localization_raw'),(AFTER,'mutated_localization_raw')]:
            values = [r.get(key) for r in raw_records]
            restored,status = restore_other(row.get(field),values)
            if status=='RESTORED_FROM_RAW': update(field,restored,status)
            elif status!='NOT_NEEDED': notes.append(field+':'+status)
            if status!='NOT_NEEDED':
                source_details.append({'id':rid,'field':field,'status':status,'members':[{'id':mid,'raw':value,'source':sources.locators.get(mid)} for mid,value in zip(members,values)]})
        if row.get(BEFORE)!=original.get(BEFORE) or row.get(AFTER)!=original.get(AFTER):
            parts = text(original.get('Mislocalization Phenotype')).split('>')
            if len(parts)==2 and location_key(parts[0])==location_key(original.get(BEFORE)) and location_key(parts[1])==location_key(original.get(AFTER)):
                update('Mislocalization Phenotype',text(row[BEFORE])+'>'+text(row[AFTER]),'Synchronize explicit before/after formatting')
            elif text(original.get('Mislocalization Phenotype')).casefold()!='no change':
                notes.append('PHENOTYPE_STRUCTURE_REQUIRES_REVIEW')
        value, variant_status, reason = normalize_variant(row['Mutation Information'],row['Mutation Type'])
        update('Mutation Information',value,reason)
        accessions = {text(r.get('uniprot_id')) for r in raw_records if ACCESSION.fullmatch(text(r.get('uniprot_id')))}
        table_acc = row.get('Accession_new')
        identity_conflict = len(accessions)>1 or (not missing(table_acc) and accessions and text(table_acc) not in accessions)
        accession = table_acc if not missing(table_acc) else next(iter(accessions)) if len(accessions)==1 else None
        if identity_conflict:
            lookup = {'status':'SOURCE_ACCESSIONS_CONFLICT','matches':[]}
        elif not missing(accession) and not ACCESSION.fullmatch(text(accession)):
            lookup = {'status':'EXISTING_ACCESSION_REQUIRES_REVIEW','matches':[]}
        else:
            lookup = client.resolve(row['Gene Name'],accession)
        residue_status = 'NOT_CHECKED'
        if lookup['status']=='UNIQUE':
            values = entry_values(lookup['matches'][0], row)
            sequence = values['Wild Sequence_new']
            if sequence:
                residue_status = sequence_check(value,sequence)
                explicit_isoform = any(not missing(r.get('isoform')) for r in raw_records)
                if explicit_isoform or residue_status in {'REFERENCE_RESIDUE_MISMATCH','POSITION_OUT_OF_RANGE'}:
                    values['Wild Sequence_new']=None
                    notes.append('SEQUENCE_WITHHELD_ISOFORM_CONTEXT' if explicit_isoform else 'SEQUENCE_WITHHELD_REFERENCE_MISMATCH')
            for field,value in values.items():
                if value is None: continue
                if missing(row.get(field)): update(field,value,'UniProt exact human identity; database annotation')
                elif text(row[field])!=text(value): notes.append(field+':EXISTING_VALUE_DIFFERS')
        audit.append({'id':rid,'PMID':row['PMID'],'mutation_status':variant_status,'mutation_note':reason,
                      'uniprot_status':lookup['status'],'candidate_accessions':[e['primaryAccession'] for e in lookup.get('matches',[])],
                      'sequence_check':residue_status,'notes':'; '.join(notes),'source_ids':members,
                      'query_url':lookup.get('response_url'),'query_time':lookup.get('fetched_at'),'database_release':lookup.get('release')})
        output.append(row)
        if i==1 or i%100==0 or i==len(rows): print(f'{i}/{len(rows)} rows; {client.calls} HTTP requests; {client.cache_hits} cache hits',flush=True)
    return output, changes, audit, source_details


def write_csv(path,headers,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as stream:
        writer = csv.writer(stream); writer.writerow(headers)
        for row in rows:
            values = [row.get(h) for h in headers]
            # CSV is a transport format; exact original strings also remain in JSON.
            values = [json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for v in values]
            values = ["'"+v if isinstance(v,str) and v.lstrip().startswith(('=','+','-','@')) else v for v in values]
            writer.writerow(values)


def parser():
    p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--sheet',help='Default: 完整记录_26列, or the unique matching sheet')
    p.add_argument('--output-dir',type=Path,default=ROOT/'data/20_enrich_legacy_table/run')
    p.add_argument('--source-audit',type=Path,help='CSV with id, PMID, 原始抽取文件, JSON定位')
    p.add_argument('--original-records',type=Path,help='Fallback Step 04 JSON with stable per-paper record order')
    p.add_argument('--dedup-map',type=Path,help='JSON containing canonical original-ID to retained-ID mapping')
    p.add_argument('--cache-dir',type=Path,default=ROOT/'data/20_enrich_legacy_table/uniprot_cache')
    mode=p.add_mutually_exclusive_group(); mode.add_argument('--online',action='store_true'); mode.add_argument('--offline',action='store_true')
    p.add_argument('--include-unreviewed',action='store_true',help='Search unreviewed entries only when no exact reviewed match exists')
    p.add_argument('--max-requests',type=int,help='HTTP request budget including retries; cached reads are free')
    p.add_argument('--limit-rows',type=int,help='Explicitly export only the first N rows, for a trial')
    p.add_argument('--timeout',type=float,default=30)
    p.add_argument('--retries',type=int,default=2)
    p.add_argument('--interval',type=float,default=0.4)
    p.add_argument('--cache-max-age-days',type=float,default=30)
    p.add_argument('--xlsx',action='store_true',help='Also export a formatted workbook using project openpyxl dependency')
    p.add_argument('--force',action='store_true',help='Replace this output run; never replaces input')
    return p


def export_workbook(out, headers, updated, changes, audit, details, *, data_only=False,
                    sheet_name='Enriched Records', filename='enriched_records.xlsx'):
    """Keep complete sequences in FASTA/JSON; split oversized Excel cells losslessly."""
    from openpyxl import Workbook, load_workbook
    long_sequences = {r['id']:r for r in updated if len(text(r.get('Wild Sequence_new')))>32767}
    def marker(rid):
        if data_only:
            return f'[Full sequence: wild_sequences.fasta / {Path(filename).with_suffix(".json").name}; id={rid}; {len(long_sequences[rid]["Wild Sequence_new"])} aa]'
        return f'[Full sequence: Long Sequences sheet / wild_sequences.fasta; id={rid}; {len(long_sequences[rid]["Wild Sequence_new"])} aa]'
    exported = [{**r, 'Wild Sequence_new':marker(r['id'])} if r['id'] in long_sequences else r for r in updated]
    exported_changes = []
    for change in changes:
        row = dict(change)
        if row['id'] in long_sequences and row['field']=='Wild Sequence_new':
            for key in ['before','after']:
                if len(text(row.get(key)))>32767: row[key]=marker(row['id'])
        exported_changes.append(row)
    with (out/'wild_sequences.fasta').open('w',encoding='utf-8',newline='\n') as stream:
        for row in updated:
            sequence = text(row.get('Wild Sequence_new'))
            if not sequence: continue
            stream.write('>'+row['id']+' accession='+text(row.get('Accession_new'))+'\n')
            for start in range(0,len(sequence),80): stream.write(sequence[start:start+80]+'\n')
    chunks = []
    for rid,row in long_sequences.items():
        sequence=row['Wild Sequence_new']; count=(len(sequence)+29999)//30000
        for i,start in enumerate(range(0,len(sequence),30000),1):
            chunks.append([rid,row['PMID'],row.get('Accession_new'),len(sequence),i,count,sequence[start:start+30000]])
    book=Workbook(write_only=True)
    N.write_sheet(book,sheet_name,headers,([r.get(h) for h in headers] for r in exported),plain=data_only)
    if not data_only:
        N.write_sheet(book,'Record Audit',list(audit[0]) if audit else ['id'],([r[h] for h in audit[0]] for r in audit))
        change_headers=['id','PMID','field','before','after','rule']
        N.write_sheet(book,'Field Changes',change_headers,([r.get(h) for h in change_headers] for r in exported_changes))
        N.write_sheet(book,'Other Sources',['id','field','status','members'],([r[h] for h in ['id','field','status','members']] for r in details))
        if chunks:
            N.write_sheet(book,'Long Sequences',['id','PMID','Accession','Sequence_Length','Chunk_Index','Chunk_Count','Sequence_Chunk'],chunks)
    book.save(out/filename)
    checked=load_workbook(out/filename,read_only=True,data_only=False)
    try:
        expected=N.data_digest([headers]+[[r.get(h) for h in headers] for r in exported])
        if N.data_digest(checked[sheet_name].values)!=expected: raise ValueError('XLSX readback mismatch')
        if data_only and checked.sheetnames != [sheet_name]: raise ValueError('Data export must have one worksheet')
        if chunks and not data_only:
            recovered=defaultdict(list)
            for r in list(checked['Long Sequences'].values)[1:]: recovered[r[0]].append((r[4],r[6]))
            if {rid:''.join(v for _,v in sorted(parts)) for rid,parts in recovered.items()} != {rid:r['Wild Sequence_new'] for rid,r in long_sequences.items()}:
                raise ValueError('Long-sequence chunk readback mismatch')
    finally: checked.close()
    return {'xlsx_readback_verified':True,'xlsx_long_sequence_records':len(long_sequences),
            'xlsx_long_sequence_ids':list(long_sequences),'complete_sequences_fasta':'wild_sequences.fasta',
            'data_only':data_only,'cell_fills':not data_only}


def main(argv=None):
    p=parser(); args=p.parse_args(argv)
    for key in ['timeout','cache_max_age_days']:
        if not 0<getattr(args,key)<float('inf'): p.error(key+' must be finite and positive')
    if not 0<=args.interval<float('inf') or args.retries<0 or (args.max_requests is not None and args.max_requests<0) or (args.limit_rows is not None and args.limit_rows<1): p.error('Invalid limits')
    src=args.input.resolve(); out=args.output_dir.resolve()
    if src.parent==out or out in src.parents: p.error('Output directory must not contain the input file')
    if out.exists() and any(out.iterdir()) and not args.force: p.error('Output exists; choose a new output directory or use --force')
    # Autodiscovery is limited to the explicitly selected workbook's adjacent audit files.
    adjacent_audit=src.parent.parent/'all_record_audit.csv'
    adjacent_dedup=src.parent/'dedup_data.json'
    audit_path=args.source_audit or (adjacent_audit if adjacent_audit.exists() else None)
    dedup_path=args.dedup_map or (adjacent_dedup if adjacent_dedup.exists() else None)
    sources=HistoricalSources(audit_path,args.original_records,dedup_path)
    input_hash=digest(src)
    headers,rows=load_rows(src,args.sheet); total=len(rows)
    if args.limit_rows: rows=rows[:args.limit_rows]
    sources.load(rows)
    generated=['enriched_records.json','field_changes.json','record_audit.json','other_source_details.json',
               'enriched_records.csv','record_audit.csv','enriched_records.xlsx','run_report.json','wild_sequences.fasta']
    protected={src,*[Path(s).resolve() for s in sources.hashes]}
    if any((out/name).resolve() in protected for name in generated):
        p.error('An output path would overwrite a historical input')
    client=UniProtClient(args.cache_dir,args.online,args.max_requests,args.timeout,args.retries,args.interval,args.cache_max_age_days,args.include_unreviewed)
    try:
        updated,changes,audit,details=process(rows,sources,client)
    finally:
        client.close()
    assert [r['id'] for r in updated]==[r['id'] for r in rows]
    if digest(src)!=input_hash: raise ValueError('Input changed while processing')
    for filename,expected in sources.hashes.items():
        if digest(Path(filename))!=expected: raise ValueError('Historical source changed: '+filename)
    out.mkdir(parents=True,exist_ok=True)
    report={'version':VERSION,'script_sha256':digest(Path(__file__)),
            'normalizer_sha256':digest(ROOT/'scripts/_internal/record_cleaning.py'),
            'input':str(src),'input_sha256':input_hash,'source_sha256':sources.hashes,
            'input_rows':total,'output_rows':len(updated),'is_partial':len(updated)!=total,
            'mode':'online' if args.online else 'offline','http_requests':client.calls,'cache_hits':client.cache_hits,
            'changed_cells':len(changes),'changed_records':len({r['id'] for r in changes}),
            'changes_by_field':dict(Counter(r['field'] for r in changes)),
            'mutation_status':dict(Counter(r['mutation_status'] for r in audit)),
            'uniprot_status':dict(Counter(r['uniprot_status'] for r in audit)),
            'other_status':dict(Counter(r['status'] for r in details)),
            'source_audit':str(audit_path) if audit_path else None,'dedup_map':str(dedup_path) if dedup_path else None,
            'final_admission':False,'reference_sequence_validation':'Only stated reference residues are checked; no renumbering',
            'created_at':datetime.now(timezone.utc).isoformat()}
    atomic_json(out/'enriched_records.json',{'columns':headers,'records':updated,'report':report})
    atomic_json(out/'field_changes.json',changes)
    atomic_json(out/'record_audit.json',audit)
    atomic_json(out/'other_source_details.json',details)
    write_csv(out/'enriched_records.csv',headers,updated)
    write_csv(out/'record_audit.csv',list(audit[0]) if audit else ['id'],audit)
    if args.xlsx:
        report.update(export_workbook(out,headers,updated,changes,audit,details))
        atomic_json(out/'enriched_records.json',{'columns':headers,'records':updated,'report':report})
    atomic_json(out/'run_report.json',report)
    print(json.dumps(report,ensure_ascii=True,indent=2))
    return 0


