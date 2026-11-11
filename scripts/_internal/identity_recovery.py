"""Recover only the identity and isoform cohorts left by step 11 (formerly 20).

Uses HGNC approved names/aliases, UniProt cross-references and explicitly named
isoforms or versioned RefSeq records. Never selects an isoform merely because
the mutation would fit. The separate reference-mismatch cohort is untouched.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('step20',ROOT/'scripts/_internal/legacy_enrichment.py')
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)
GREEK={'伪':'alpha','尾':'beta','纬':'gamma','未':'delta','蔚':'epsilon','胃':'theta','魏':'kappa','位':'lambda','渭':'mu','碌':'mu','蟽':'sigma','蟿':'tau','味':'zeta'}


def norm(value):
    value=unicodedata.normalize('NFKC',E.text(value)).casefold()
    for a,b in GREEK.items():value=value.replace(a,b)
    return re.sub(r'[^a-z0-9]','',value)


def components(value):
    """Split explicit name aliases, not arbitrary words or hyphenated constructs."""
    value=E.text(value)
    return list(dict.fromkeys([value]+[p.strip() for p in re.split(r'[();/]',value) if p.strip()]))


def source_aliases(work):
    """Read explicit local alias definitions, never infer identity from residues."""
    if not work.get('source_path'):
        return []
    path=(ROOT/work['source_path']).resolve()
    if not path.is_relative_to(ROOT) or E.digest(path)!=work.get('source_sha256'):
        raise ValueError('Identity source path or SHA256 mismatch')
    source=path.read_text(encoding='utf-8-sig')
    source=re.split(r'(?im)^\s*#{1,6}\s*(?:references|bibliography)\s*$',source)[0]
    name=E.text(work['row']['Gene Name'])
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9-]{1,25}',name):
        return []
    token=r'([A-Za-z][A-Za-z0-9-]{1,25})'
    patterns=[r'(?<![A-Za-z0-9])'+re.escape(name)+r'\s*\(\s*'+token+r'\s*\)']
    # A lower-case h is expanded only when the article defines it as human.
    if name.startswith('h'):
        patterns += [r'(?<![A-Za-z0-9])'+re.escape(name)+r'\s*,\s*human\s+'+token+r'\b',
                     r'\bhuman\s+'+token+r'\s*\(\s*'+re.escape(name)+r'\s*\)']
    return [{'alias':m[1], 'quote':m[0], 'source_path':work['source_path'],
             'source_sha256':work['source_sha256'], 'basis':'explicit source alias definition'}
            for pattern in patterns for m in re.finditer(pattern,source,re.I)]


def protein_names(entry):
    """Protein product names only: gene names cannot distinguish shared loci."""
    description=entry.get('proteinDescription',{})
    names=[]
    for item in [description.get('recommendedName',{}),*description.get('alternativeNames',[])]:
        names.append(item.get('fullName',{}).get('value',''))
        names.extend(v.get('value','') for v in item.get('shortNames',[]))
    return {norm(n) for n in names if n}


class References:
    def __init__(self, hgnc, old_cache, new_cache, online, max_requests=None, max_age_days=30, include_unreviewed=False):
        self.max_age_days=max_age_days; self.include_unreviewed=include_unreviewed
        self.hgnc={d['hgnc_id']:d for d in E.read_json(hgnc)['response']['docs'] if d.get('status')=='Approved'}
        self.primary=defaultdict(set);self.alias=defaultdict(set);self.entries={};self.entry_sources={}
        for hid,doc in self.hgnc.items():
            self.primary[norm(doc['symbol'])].add(hid)
            for field in ['symbol','name','alias_symbol','alias_name','prev_symbol','prev_name']:
                values=doc.get(field) or [];values=values if isinstance(values,list) else [values]
                for value in values:self.alias[norm(value)].add(hid)
        self.old_cache=old_cache;self.new_cache=new_cache;self.online=online;self.max_requests=max_requests
        self.client=E.UniProtClient(new_cache/'uniprot_json',online=online,interval=0.5,max_requests=max_requests,max_age_days=max_age_days,include_unreviewed=include_unreviewed)
        self.session=requests.Session();self.text_requests=0;self.last_text=0
        for directory in [old_cache,new_cache/'uniprot_json']:
            for path in sorted(directory.glob('*.json')):
                envelope=E.read_json(path)
                if envelope.get('status')!='OK' or envelope.get('schema')!='uniprot-request-cache-v1' or not self.fresh(envelope):continue
                request=envelope.get('request',{})
                if not request.get('url','').startswith('https://rest.uniprot.org/uniprotkb/'):continue
                key=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
                if path.stem!=key:continue
                payload=envelope.get('payload') or {}
                for entry in payload.get('results',[payload]):
                    if entry.get('primaryAccession') and entry.get('organism',{}).get('taxonId')==9606:
                        self.add_entry(entry,envelope.get('response_url'),envelope.get('release'),envelope.get('fetched_at'))

    def fresh(self,envelope):
        try:
            age=(datetime.now(timezone.utc)-datetime.fromisoformat(envelope["fetched_at"])).total_seconds()/86400
            return 0 <= age <= self.max_age_days
        except (KeyError,TypeError,ValueError):
            return False

    def entry_hgnc(self,entry):
        ids={x['id'] for x in entry.get('uniProtKBCrossReferences',[]) if x.get('database')=='HGNC' and x.get('id') in self.hgnc}
        if not ids:
            for gene in entry.get('genes',[]):ids.update(self.primary.get(norm(gene.get('geneName',{}).get('value')),set()))
        return ids

    def add_entry(self,entry,url,release,fetched_at=None):
        accession=entry['primaryAccession'];self.entries[accession]=entry
        self.entry_sources[accession]={'url':url,'release':release,'fetched_at':fetched_at}
        ids=self.entry_hgnc(entry)
        for value in E.N.HumanUniProtLookup.exact_names(entry):self.alias[norm(value)].update(ids)

    def get_entry(self,accession):
        accession=accession.split('-')[0]
        if accession in self.entries:return self.entries[accession],None
        for entry in list(self.entries.values()):
            if accession in entry.get('secondaryAccessions',[]):return entry,None
        if self.max_requests is not None:self.client.max_requests=max(0,self.max_requests-self.text_requests)
        response=self.client.request(accession+'.json')
        if response['status']!='OK':return None,response['status']
        entry=response['payload']
        if entry.get('organism',{}).get('taxonId')!=9606:return None,'NOT_HUMAN'
        if entry.get('primaryAccession')!=accession and accession not in entry.get('secondaryAccessions',[]):return None,'ACCESSION_REDIRECT_CONFLICT'
        self.add_entry(entry,response.get('response_url'),response.get('release'),response.get('fetched_at'))
        self.entries[accession]=entry
        return entry,None

    def name_ids(self,value):
        key=norm(value)
        return self.primary.get(key) or self.alias.get(key,set())

    def identify(self,work):
        evidence=[];sets=[]
        gene=work['row']['Gene Name']; originals=[o['record'] for o in work['originals']]
        for value in components(gene):
            ids=set(self.name_ids(value))
            if ids:sets.append(ids);evidence.append({'value':value,'hgnc_ids':sorted(ids),'basis':'HGNC/UniProt name'})
        for proof in source_aliases(work):
            ids=set(self.name_ids(proof['alias']))
            if ids:
                sets.append(ids);evidence.append({**proof,'hgnc_ids':sorted(ids)})
        for record in originals:
            h=E.text(record.get('hgnc_id'))
            if h in self.hgnc:sets.append({h});evidence.append({'value':h,'hgnc_ids':[h],'basis':'source HGNC ID'})
        # Observed-protein names can corroborate/disambiguate an existing name.
        # Alone they are used only where the original gene label occurs explicitly
        # as an alias in that same observed-protein description.
        observed_sets=[]
        for record in originals:
            observed=record.get('observed_protein')
            linked=norm(gene) in {norm(p) for p in components(observed)}
            for value in components(observed):
                ids=set(self.name_ids(value))
                if ids and (sets or linked):
                    observed_sets.append(ids);evidence.append({'value':value,'hgnc_ids':sorted(ids),'basis':'source observed-protein corroboration'})
        if sets:
            candidates=set.intersection(*sets)
            if not candidates:return None,'GENE_LABELS_CONFLICT',evidence
            for ids in observed_sets:
                overlap=candidates&ids
                if not overlap:return None,'OBSERVED_IDENTITY_CONFLICT',evidence
                candidates=overlap
        elif observed_sets:candidates=set.intersection(*observed_sets)
        else:return None,'NO_EXACT_AUTHORITY_NAME',evidence
        if len(candidates)!=1:return None,'NAME_STILL_AMBIGUOUS',evidence
        hid=next(iter(candidates));doc=self.hgnc[hid]
        declared={E.text(r.get('uniprot_id')).split('-')[0] for r in originals if E.ACCESSION.fullmatch(E.text(r.get('uniprot_id')))}
        existing=E.text(work['row'].get('Accession_new'))
        if not E.missing(existing):declared.add(existing.split('-')[0])
        accessions=set(doc.get('uniprot_ids') or [])
        if declared:
            resolved=set()
            for accession in sorted(declared):
                entry,error=self.get_entry(accession)
                if error:return None,'SOURCE_ACCESSION_UNRESOLVED',evidence
                resolved.add(entry['primaryAccession'])
                evidence.append({'declared_accession':accession,'primary_accession':entry['primaryAccession'],'basis':'UniProt accession mapping'})
            if len(resolved)!=1:return None,'SOURCE_ACCESSIONS_CONFLICT',evidence
            canonical=next(iter(resolved))
            if canonical not in accessions:return None,'SOURCE_ACCESSION_CONFLICT',evidence
            accessions={canonical}
        entries=[]
        for accession in sorted(accessions):
            entry,error=self.get_entry(accession)
            if error:return None,'ENTRY_'+error,evidence
            if hid not in self.entry_hgnc(entry):return None,'HGNC_UNIPROT_CONFLICT',evidence
            entries.append(entry)
        reviewed=[e for e in entries if e.get('entryType','').startswith('UniProtKB reviewed')]
        entries=(reviewed or entries) if self.include_unreviewed else reviewed
        if not entries:return None,'UNREVIEWED_NOT_ALLOWED',evidence
        unique={e['primaryAccession']:e for e in entries}
        if len(unique)>1:
            # Superscript state labels do not change protein product identity.
            # Each merged member must independently select the same product.
            member_matches=[]
            for original in originals:
                observed=re.sub(r'\^\{?[A-Za-z]+\}?','',E.text(original.get('observed_protein')))
                names={norm(v) for v in components(observed)}
                member_matches.append({acc for acc,e in unique.items() if names & protein_names(e)})
            selected=set.intersection(*member_matches) if member_matches else set()
            if len(selected)==1:
                accession=next(iter(selected))
                evidence.append({'basis':'observed protein product name, distinct from gene symbol',
                                 'observed_proteins':[o.get('observed_protein') for o in originals],
                                 'candidate_accessions':sorted(unique),'selected_accession':accession})
                unique={accession:unique[accession]}
        if len(unique)!=1:return None,'MULTIPLE_OR_MISSING_AUTHORITY_ENTRIES',evidence
        entry=next(iter(unique.values()))
        evidence.append({'hgnc_id':hid,'symbol':doc['symbol'],'uniprot':entry['primaryAccession'],'basis':'HGNC UniProt cross-reference verified in human entry'})
        return entry,'IDENTITY_CONFIRMED',evidence

    def fetch_text(self,url,params=None):
        identity={'url':url,'params':params or {}}
        key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
        path=self.new_cache/'text'/f'{key}.json'
        if path.exists():
            cached=E.read_json(path)
            if cached.get('request')==identity and cached.get('status')=='OK' and self.fresh(cached):return cached
            if not self.online:return {'status':'STALE_CACHE','request':identity}
        if not self.online:return {'status':'OFFLINE_CACHE_MISS','request':identity}
        for attempt in range(3):
            if self.max_requests is not None and self.client.calls+self.text_requests>=self.max_requests:
                return {'status':'REQUEST_LIMIT','request':identity}
            time.sleep(max(0,0.6-(time.monotonic()-self.last_text)))
            self.last_text=time.monotonic();self.text_requests+=1
            try:
                response=self.session.get(url,params=params,timeout=(15,60))
                response.raise_for_status()
                result={'status':'OK','text':response.text,'url':response.url,'request':identity,
                        'fetched_at':datetime.now(timezone.utc).isoformat(),
                        'release':response.headers.get('X-UniProt-Release')}
                E.atomic_json(path,result);return result
            except requests.RequestException as exc:
                if attempt==2:return {'status':'REQUEST_FAILED','request':identity,'error':str(exc)}
                time.sleep(2**attempt)

    def isoforms(self,entry):
        return [iso for c in entry.get('comments',[]) if c.get('commentType')=='ALTERNATIVE PRODUCTS' for iso in c.get('isoforms',[])]

    def iso_sequence(self,entry,isoid):
        iso=next((i for i in self.isoforms(entry) if isoid in i.get('isoformIds',[])),None)
        if iso is None:return None,'ISOFORM_NOT_IN_ENTRY',{}
        if iso.get('isoformSequenceStatus')=='Displayed':
            return entry.get('sequence',{}).get('value'),'UNIPROT_DISPLAYED_ISOFORM',{'sequence_reference':isoid,'url':self.entry_sources.get(entry['primaryAccession'],{}).get('url')}
        if iso.get('isoformSequenceStatus')!='Described':return None,'ISOFORM_SEQUENCE_UNAVAILABLE',{}
        result=self.fetch_text(f'https://rest.uniprot.org/uniprotkb/{isoid}.fasta')
        if result['status']!='OK':return None,result['status'],result
        lines=result['text'].splitlines()
        if not lines or len([x for x in lines if x.startswith('>')])!=1:return None,'INVALID_FASTA',{}
        bits=lines[0].split('|')
        if len(bits)<3 or bits[1]!=isoid:return None,'FASTA_ID_MISMATCH',{}
        sequence=''.join(lines[1:]).strip()
        if not re.fullmatch('[A-Z]+',sequence):return None,'INVALID_FASTA_SEQUENCE',{}
        return sequence,'UNIPROT_NAMED_ISOFORM',{'sequence_reference':isoid,'url':result['url'],'retrieved_at':result['fetched_at']}

    def refseq_sequence(self,reference,entry):
        db='nuccore' if reference.startswith('NM_') else 'protein'
        result=self.fetch_text('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',{'db':db,'id':reference,'rettype':'gb','retmode':'xml'})
        if result['status']!='OK':return None,result['status'],result
        try:items=ET.fromstring(result['text']).findall('.//GBSeq')
        except ET.ParseError:return None,'INVALID_REFSEQ_XML',{}
        if len(items)!=1:return None,'NONUNIQUE_REFSEQ_RESPONSE',{}
        item=items[0];version=item.findtext('GBSeq_accession-version','')
        if ('.' in reference and version!=reference) or version.split('.')[0]!=reference.split('.')[0]:return None,'REFSEQ_VERSION_MISMATCH',{}
        if item.findtext('GBSeq_organism')!='Homo sapiens':return None,'REFSEQ_NOT_HUMAN',{}
        allowed=set()
        for hid in self.entry_hgnc(entry):
            doc=self.hgnc[hid];allowed.add(norm(doc['symbol']))
            for key in ['alias_symbol','prev_symbol']:allowed.update(norm(v) for v in doc.get(key) or [])
        xrefs={v.split('.')[0] for x in entry.get('uniProtKBCrossReferences',[]) if x.get('database')=='RefSeq'
               for v in [x['id']]+[p['value'] for p in x.get('properties',[]) if p['key']=='NucleotideSequenceId']}
        cds=[];gene_ids=set()
        for feature in item.findall('.//GBFeature'):
            qualifiers=defaultdict(list)
            for q in feature.findall('.//GBQualifier'):qualifiers[q.findtext('GBQualifier_name')].append(q.findtext('GBQualifier_value',''))
            gene_ids.update(norm(v) for v in qualifiers['gene'])
            if feature.findtext('GBFeature_key')=='CDS':cds.append(qualifiers)
        if gene_ids and not gene_ids&allowed:return None,'REFSEQ_GENE_CONFLICT',{}
        if not gene_ids&allowed and reference.split('.')[0] not in xrefs:return None,'REFSEQ_GENE_UNBOUND',{}
        if db=='nuccore':
            translated=[q for q in cds if q['translation'] and q['protein_id']]
            if len(translated)!=1:return None,'REFSEQ_CDS_NOT_UNIQUE',{}
            sequence=translated[0]['translation'][0].replace(' ','').replace('\n','')
            protein=translated[0]['protein_id'][0]
        else:sequence=item.findtext('GBSeq_sequence','').upper();protein=version
        if not re.fullmatch('[A-Z]+',sequence):return None,'REFSEQ_SEQUENCE_INVALID',{}
        return sequence,'EXPLICIT_REFSEQ_VERSION',{'sequence_reference':protein,'transcript_reference':version if db=='nuccore' else None,'url':result['url'],'retrieved_at':result['fetched_at']}

    def choose_one_isoform(self,label,entry):
        """Use explicit IDs or labels; never rank versions by residue compatibility."""
        label=E.text(label)
        ids=re.findall(r'\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})-[1-9]\d*\b',label)
        if ids:
            if len(set(ids))!=1:return None,'MULTIPLE_ISOFORM_REFERENCES',{}
            return self.iso_sequence(entry,ids[0])
        standardized=re.sub(r'\b(N[MP])[_ ]?(\d+)',r'\1_\2',label)
        refs=re.findall(r'\bN[MP]_\d+(?:\.\d+)?\b',standardized)
        if refs:
            if len(set(refs))!=1:return None,'MULTIPLE_REFSEQ_REFERENCES',{}
            return self.refseq_sequence(refs[0],entry)
        ensembl=re.findall(r'\bENS[TP]\d+(?:\.\d+)?\b',label)
        if ensembl:
            hits=set()
            for x in entry.get('uniProtKBCrossReferences',[]):
                if x.get('database')!='Ensembl':continue
                values=[x['id']]+[p['value'] for p in x.get('properties',[]) if p['key']=='ProteinId']
                if any(v==r or ('.' not in r and v.split('.')[0]==r) for v in values for r in ensembl) and x.get('isoformId'):hits.add(x['isoformId'])
            if len(hits)==1:return self.iso_sequence(entry,next(iter(hits)))
            return None,'ENSEMBL_VERSION_OR_ISOFORM_UNRESOLVED',{}
        hits=set()
        # Reject explicit domain-only/engineered fragments; they are not natural isoforms.
        if re.search(r'\bfragment\b|repeat domain|repeat region|\bdel\d+\b',label,re.I):return None,'FRAGMENT_OR_CONSTRUCT_CONTEXT',{}
        parts=[norm(v) for v in components(label)]
        for iso in self.isoforms(entry):
            aliases=[iso.get('name',{}).get('value','')]+[s['value'] for s in iso.get('synonyms',[])]
            for alias in aliases:
                key=norm(alias)
                if not key:continue
                if key in parts or norm('isoform '+alias) in parts:
                    hits.update(iso.get('isoformIds',[]))
                elif len(key)>=3 and re.search(r'(?<![a-z0-9])'+re.escape(alias.casefold())+r'(?![a-z0-9])',label.casefold()):
                    hits.update(iso.get('isoformIds',[]))
        if len(hits)==1:return self.iso_sequence(entry,next(iter(hits)))
        if len(hits)>1:return None,'MULTIPLE_ISOFORMS_IN_LABEL',{}
        # A label that only repeats the gene is unspecified, as in the original
        # step-20 records with an empty isoform field; retain canonical-reference status.
        names={norm(v) for v in E.N.HumanUniProtLookup.exact_names(entry)}
        generic=re.sub(r'^(?:full[- ]length\s+)?(?:human\s+)?','',label,flags=re.I)
        if norm(generic) in names:return entry.get('sequence',{}).get('value'),'GENERIC_NAME_CANONICAL_REFERENCE',{'sequence_reference':entry['primaryAccession'],'url':self.entry_sources.get(entry['primaryAccession'],{}).get('url')}
        return None,'ISOFORM_LABEL_UNRESOLVED',{}

    def choose_sequence(self,work,entry):
        per_member=[]
        for original in work['originals']:
            raw=original['record'];label=E.text(raw.get('isoform'))
            explicit=E.text(raw.get('uniprot_id'))
            if E.ACCESSION.fullmatch(explicit) and '-' in explicit:
                if not E.missing(label) and re.search(r'\b(?:NM_|NP_|ENS[TP])',label):
                    # Conflicting independent references are not resolved by priority.
                    a=self.choose_one_isoform(explicit,entry);b=self.choose_one_isoform(label,entry)
                    if not a[0] or not b[0] or a[0]!=b[0]:return None,'SOURCE_REFERENCES_DISAGREE',{}
                    per_member.append(a);continue
                label=explicit
            if E.missing(label):
                per_member.append((entry.get('sequence',{}).get('value'),'CANONICAL_REFERENCE',{'sequence_reference':entry['primaryAccession'],'url':self.entry_sources.get(entry['primaryAccession'],{}).get('url')}))
            else:per_member.append(self.choose_one_isoform(label,entry))
        if any(not item[0] for item in per_member):
            failed=next(item for item in per_member if not item[0]);return failed
        if len({item[0] for item in per_member})!=1:return None,'MERGED_SOURCE_SEQUENCE_REFERENCES_DIFFER',{}
        # Even equal sequences retain all independently cited source references.
        result=per_member[0]
        return result[0],result[1],{**result[2],'member_references':[i[2] for i in per_member]}

    def close(self):self.client.close();self.session.close()


def run(args):
    prior=E.read_json(args.input_dir/'enriched_records.json')
    audit=E.read_json(args.input_dir/'record_audit.json');audit_by={a['id']:a for a in audit}
    workset=E.read_json(args.workset)
    if E.digest(args.input_dir/'enriched_records.json')!=workset['input_sha256']:raise ValueError('Workset is stale')
    if args.output_dir.resolve()==args.input_dir.resolve():raise ValueError('Output cannot replace the input directory')
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.force:raise ValueError('Use a new empty output directory or --force for this output run')
    references=References(args.hgnc,args.old_cache,args.cache_dir,args.online,getattr(args,'max_requests',None),getattr(args,'cache_max_age_days',30),getattr(args,'include_unreviewed',False))
    rows=copy.deepcopy(prior['records']);by_id={r['id']:r for r in rows}
    recovery=[];changes=[]
    try:
        for index,work in enumerate(workset['records'],1):
            row=by_id[work['row']['id']];rid=row['id'];rec={'id':rid,'PMID':row['PMID'],'cohort':work['cohort']}
            entry,status,evidence=references.identify(work)
            rec.update(identity_status=status,identity_evidence=evidence,sequence_status='IDENTITY_PENDING',sequence_reference=None,reference_check='NOT_CHECKED')
            if entry:
                values=E.entry_values(entry,row);seq,seq_status,proof=references.choose_sequence(work,entry)
                rec.update(sequence_status=seq_status,sequence_reference=proof.get('sequence_reference'),sequence_evidence=proof)
                if seq:
                    check=E.sequence_check(row['Mutation Information'],seq);rec['reference_check']=check
                    if check in {'REFERENCE_RESIDUE_MISMATCH','POSITION_OUT_OF_RANGE'}:
                        seq=None;rec['sequence_status']='REFERENCE_MISMATCH_REMAINS'
                values['Wild Sequence_new']=seq
                # Only identity cohort can acquire new metadata; isoform cohort
                # only fills its sequence and preserves existing database columns.
                allowed=E.ENRICH_FIELDS if work['cohort']=='identity' else ['Wild Sequence_new']
                for field in allowed:
                    value=values.get(field)
                    if value is not None and E.missing(row.get(field)):
                        changes.append({'id':rid,'PMID':row['PMID'],'field':field,'before':row.get(field),'after':value,'rule':'Step21 '+status+'; '+rec['sequence_status']})
                        row[field]=value
                rec['matched_accession']=entry['primaryAccession']
                a=audit_by[rid]
                a['uniprot_status']='UNIQUE';a['candidate_accessions']=[entry['primaryAccession']]
                a['query_url']=references.entry_sources.get(entry['primaryAccession'],{}).get('url')
                a['database_release']=references.entry_sources.get(entry['primaryAccession'],{}).get('release')
                a['query_time']=references.entry_sources.get(entry['primaryAccession'],{}).get('fetched_at')
                if seq:
                    a['sequence_check']=rec['reference_check']
                    a['notes']='; '.join(v for v in a['notes'].split('; ') if not v.startswith('SEQUENCE_WITHHELD'))
                a['notes']+=('; ' if a['notes'] else '')+'STEP21_'+rec['sequence_status']
            rec['filled_fields']=[c['field'] for c in changes if c['id']==rid]
            recovery.append(rec)
            if index==1 or index%50==0 or index==len(workset['records']):
                print(f'{index}/{len(workset["records"])} targets; added {sum(c["field"]=="Wild Sequence_new" for c in changes)} sequences; {references.client.calls+references.text_requests} HTTP requests',flush=True)
        assert [r['id'] for r in rows]==[r['id'] for r in prior['records']]
        targeted={w['row']['id'] for w in workset['records']}
        assert all(new==old for new,old in zip(rows,prior['records']) if old['id'] not in targeted)
        out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
        report={'version':'1.0.0','script_sha256':E.digest(Path(__file__)),'input_sha256':workset['input_sha256'],'hgnc_sha256':E.digest(args.hgnc),
                'cohorts':dict(Counter(w['cohort'] for w in workset['records'])),
                'input_rows':len(rows),'output_rows':len(rows),'outside_cohorts_unchanged':True,
                'added_by_field':dict(Counter(c['field'] for c in changes)),
                'sequence_recovery_by_cohort':dict(Counter(r['cohort'] for r in recovery if 'Wild Sequence_new' in r['filled_fields'])),
                'identity_status':dict(Counter(r['identity_status'] for r in recovery)),
                'sequence_status':dict(Counter(r['sequence_status'] for r in recovery)),
                'http_requests':references.client.calls+references.text_requests,
                'coverage':{field:sum(not E.missing(r.get(field)) for r in rows) for field in E.ENRICH_FIELDS},
                'created_at':datetime.now(timezone.utc).isoformat(),'final_admission':False}
        E.atomic_json(out/'recovery_audit.json',recovery)
        E.atomic_json(out/'field_changes.json',changes)
        E.atomic_json(out/'record_audit.json',audit)
        details=E.read_json(args.input_dir/'other_source_details.json')
        E.atomic_json(out/'other_source_details.json',details)
        E.atomic_json(out/'enriched_records.json',{'columns':prior['columns'],'records':rows,'report':report})
        E.write_csv(out/'enriched_records.csv',prior['columns'],rows)
        recovery_headers=list(dict.fromkeys(k for r in recovery for k in r))
        E.write_csv(out/'recovery_audit.csv',recovery_headers,recovery)
        if args.xlsx:
            report.update(E.export_workbook(out,prior['columns'],rows,changes,audit,details))
            from openpyxl import load_workbook
            book=load_workbook(out/'enriched_records.xlsx')
            E.N.write_sheet(book,'Identity Isoform Recovery',recovery_headers,([r.get(k) for k in recovery_headers] for r in recovery))
            book.save(out/'enriched_records.xlsx');book.close()
        E.atomic_json(out/'run_report.json',report)
        E.atomic_json(out/'enriched_records.json',{'columns':prior['columns'],'records':rows,'report':report})
        print(json.dumps(report,ensure_ascii=True,indent=2))
    finally:references.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=ROOT/'data/20_enrich_legacy_table/online_full_20261109')
    p.add_argument('--workset',type=Path,default=ROOT/'data/21_identity_isoform_recovery/workset.json')
    p.add_argument('--hgnc',type=Path,default=ROOT/'data/21_identity_isoform_recovery/reference/hgnc_complete_set.json')
    p.add_argument('--old-cache',type=Path,default=ROOT/'data/20_enrich_legacy_table/uniprot_cache')
    p.add_argument('--cache-dir',type=Path,default=ROOT/'data/21_identity_isoform_recovery/reference/cache')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--online',action='store_true');p.add_argument('--xlsx',action='store_true')
    p.add_argument('--force',action='store_true',help='Replace this output run, never the input directory')
    args=p.parse_args();run(args)


