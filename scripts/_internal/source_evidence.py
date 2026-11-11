"""Shared source paragraphs, bounded JSON parsing and quotation typography.

Kept separate from command-line stages so cached source IDs and response parsing
remain identical after retiring the historical admission pipeline.
"""
import copy
import html
import json
import re
import unicodedata

APPENDIX_HEADING = re.compile(r'(?im)^(?:#{1,6}\s*)?(?:(?:figure|fig\.?|table)\s+(?:s?\d+\b|legends\b)|supplement(?:ary|al)\s+(?:methods|results)\b)[^\n]*$')


def paragraphs(text):
    cutoff = re.search(r'(?im)^\s*#{1,6}\s*(?:references|bibliography)\s*$', text)
    body = text[:cutoff.start()] if cutoff else text
    result = {}
    for match in re.finditer(r'\S[\s\S]*?(?=\n\s*\n|\Z)', body):
        value = match[0].strip()
        if not value or value.startswith('[IMAGE_REMOVED'): continue
        for offset in range(0, len(value), 6000):
            identifier = f'P{len(result)+1:04d}'
            result[identifier] = {'text': value[offset:offset+6000],
                                  'source_offset': match.start()+offset}
    return result


def complete_text_units(text,reviewer):
    """Preserve every legacy P identifier; append post-reference captions as A IDs."""
    legacy=reviewer.paragraphs(text);units=copy.deepcopy(legacy)
    cutoff=re.search(r'(?im)^\s*#{1,6}\s*(?:references|bibliography)\s*$',text)
    start=None
    if cutoff:
        headings=list(APPENDIX_HEADING.finditer(text,cutoff.end()))
        if headings:start=headings[0].start()
    if start is not None:
        index=0
        for match in re.finditer(r'\S[\s\S]*?(?=\n\s*\n|\Z)',text[start:]):
            value=match[0].strip()
            if not value or value.startswith('[IMAGE_REMOVED'):continue
            for offset in range(0,len(value),6000):
                index+=1;units[f'A{index:04d}']={'text':value[offset:offset+6000],
                    'source_offset':start+match.start()+offset,'section':'POST_REFERENCE_CAPTIONS_OR_SUPPLEMENTS'}
    assert {k:units[k] for k in legacy}==legacy
    return units,{'legacy_units':len(legacy),'appended_units':len(units)-len(legacy),
        'reference_cutoff':cutoff.start() if cutoff else None,'appendix_start':start,
        'source_characters':len(text),'empty_source':not bool(text.strip())}


def normalize_quote_text(value):
    """Normalize typography, not meaning, spelling, numbers or negation."""
    value = unicodedata.normalize('NFC', html.unescape(value))
    value = value.translate(str.maketrans({'\u2018':"'",'\u2019':"'",'\u201c':'"','\u201d':'"',
        '\u2010':'-','\u2011':'-','\u2013':'-','\u2014':'-','\ufb01':'fi','\ufb02':'fl'}))
    value = re.sub(r'\\(?:mathrm|mathbf|text|textrm)\{([^{}]*)\}',r'\1',value)
    value = value.replace('$','')
    value = re.sub(r'\*\*([^*]+)\*\*',r'\1',value)
    value = re.sub(r'`([^`]+)`',r'\1',value)
    return re.sub(r'\s+',' ',value).strip()


def parse_review_response(content, expected_count):
    """Repair only omitted object closers at record/list boundaries; log repairs.

    Never alter a string, value, verdict or citation. Other syntax errors fail.
    """
    repairs = []
    for _ in range(expected_count * 4 + 4):
        try:
            payload = json.loads(content)
            break
        except json.JSONDecodeError as exc:
            if exc.msg == 'Extra data' and re.fullmatch(r'[}\s]+',content[exc.pos:]):
                repairs.append({'operation':'remove_trailing_object_closers','offset':exc.pos})
                content = content[:exc.pos]
                continue
            if (exc.msg == "Expecting ',' delimiter" and content[exc.pos:exc.pos+1] == '}'
                    and re.match(r'}\s*(?:,\s*{"record_id"|\])',content[exc.pos:])):
                content = content[:exc.pos] + content[exc.pos+1:]
                repairs.append({'operation':'remove_unmatched_object_closer','offset':exc.pos})
                continue
            if (exc.msg == 'Expecting property name enclosed in double quotes' and
                    content.startswith('{"record_id"', exc.pos) and content[exc.pos-1:exc.pos] == ','):
                position = exc.pos-1
            elif exc.msg == "Expecting ',' delimiter" and content[exc.pos:exc.pos+1] == ']':
                position = exc.pos
            else:
                raise
            content = content[:position] + '}' + content[position:]
            repairs.append({'operation':'insert_missing_object_closer','offset':position})
    else:
        raise ValueError('JSON structure repair limit exceeded')
    for review in payload['reviews']:
        criteria = review.get('criteria')
        if isinstance(criteria,dict):
            for key in ['requires_figure_inspection','reason','corrected_fields','field_evidence',
                        'admission','admission_evidence']:
                if key not in review and key in criteria:
                    review[key] = criteria.pop(key)
                    repairs.append({'operation':'lift_misnested_review_field','record_id':review.get('record_id'),'field':key})
    return payload, repairs
