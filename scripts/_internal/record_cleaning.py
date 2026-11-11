"""Step 05: JSON aggregation with the frozen 1.1.0 normalization engine."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import chain
from urllib.parse import urlencode
from urllib.request import Request, urlopen

VERSION = "1.1.0"
BASE_COLUMNS = [
    "id", "PMID", "Gene Name", "Mutation Type", "Mutation Information",
    "Mislocalization Phenotype", "Mislocalization Phenotype\uff08before\uff09",
    "Mislocalization Phenotype\uff08after\uff09", "Subcellular_Location_Values",
    "Accession_new", "Reviewed", "Cell Line", "Disease Information", "Title",
    "Variant Clinical Classification", "Wild Sequence_new", "Experimental_Method",
    "Evidence_Type", "Tag_Strategy", "Expression_Regime", "Quantitative_Reported",
    "Redistribution_Extent", "No_Change_Basis", "Outcome_Provenance", "Cell_System_Type",
]
DETAIL = "Localization Change Detail"
BLANK_COLUMNS = {"Subcellular_Location_Values", "Accession_new", "Reviewed",
                 "Variant Clinical Classification", "Wild Sequence_new"}
REQUIRED = {"PMID", "Mutation Information", "Mislocalization Phenotype"}
UNKNOWN = {"", "unknown", "not_reported", "undetermined", "null", "n/a", "na"}
AA = dict(zip("ARNDCEQGHILKMFPSTWYV", ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu",
              "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr",
              "Trp", "Tyr", "Val"]))
A = "(?:" + "|".join(AA.values()) + ")"
POS, SEQ = A + r"[1-9]\d*", "(?:" + A + ")+"
BODY = (r"(?:" + POS + "(?:_" + POS + ")?(?:delins" + SEQ + "|ins" + SEQ + "|del|dup)|"
        + POS + "(?:" + A + r")?fs(?:(?:Ter|\*)(?:[1-9]\d*|\?))?|"
        + POS + r"(?:Ter|\*|" + A + r"|=))")
P_RE = re.compile(r"(?<![A-Za-z0-9])p\.(?:\(" + BODY + r"\)|" + BODY + r")(?![A-Za-z0-9_?])")
TYPE_ALIASES = {
    "missense": "Missense", "nonsense": "Nonsense", "frameshift": "Frameshift",
    "deletion": "Deletion", "insertion": "Insertion", "duplication": "Duplication",
    "internal tandem duplication": "Duplication", "delins": "Delins",
    "deletion with insertion": "Delins", "splice-site": "Splice-site",
    "splice site": "Splice-site", "synonymous": "Synonymous", "expansion": "Repeat expansion",
    "polyglutamine expansion": "Repeat expansion", "repeat expansion": "Repeat expansion",
    "multiple": "Complex", "complex": "Complex", "other": "Other", "unresolved": "Unresolved",
    "non-sequence perturbation": "Non-sequence perturbation",
}
TYPES = set(TYPE_ALIASES.values()) | {"Wild-type"}

# \u63a8\u8350\u540d\u79f0\u3001CVCL\u7f16\u53f7\u3001\u4fdd\u5b88\u540c\u4e49\u8bcd\u767d\u540d\u5355\u3002\u5e76\u975eCellosaurus\u5168\u5e93\uff0c\u4e0d\u63a8\u65ad\u672a\u77e5\u540d\u79f0\u3002
CELL_LINES = [
    ("HEK293", "CVCL_0045", ["HEK-293", "HEK 293", "HEK/293", "(HEK)293", "HEK;293", "293", "293 HEK", "Graham 293", "Human Embryonic Kidney 293"]),
    ("HEK293T", "CVCL_0063", ["HEK-293T", "HEK 293T", "HEK-293-T", "HEK 293 T", "293T", "293-T"]),
    ("HeLa", "CVCL_0030", ["He La", "He-La", "HeLa-CCL2", "Henrietta Lacks cells"]),
    ("HeLa Kyoto", "CVCL_1922", ["HeLa-Kyoto", "HeLa (Kyoto)", "HeLa(Kyoto)", "HeLa 'Kyoto'"]),
    ("COS-7", "CVCL_0224", ["COS7", "COS 7"]),
    ("COS-1", "CVCL_0223", ["COS1", "COS 1"]),
    ("CHO-K1", "CVCL_0214", ["CHO K1", "CHOK1"]),
    ("CHO", "CVCL_0213", ["Chinese Hamster Ovary", "CHO-ori"]),
    ("U2OS", "CVCL_0042", ["U-2 OS", "U-2OS", "U-2-OS", "U2-OS"]),
    ("tsA201", "CVCL_2737", ["tsA-201"]),
    ("NIH 3T3", "CVCL_0594", ["NIH/3T3", "NIH-3T3", "NIH3T3", "3T3NIH"]),
    ("SH-SY5Y", "CVCL_0019", ["SHSY5Y", "SHSY-5Y", "SK-SH-SY5Y"]),
    ("Hep-G2", "CVCL_0027", ["HepG2", "Hep G2"]),
    ("MDCK", "CVCL_0422", ["MDCK (NBL-2)", "MDCK(NBL-2)", "NBL-2"]),
    ("MDCK-II", "CVCL_0424", ["MDCK II", "MDCKII", "MDCK2", "MDCK-2", "MDCK Type II"]),
    ("BHK-21", "CVCL_1914", ["BHK 21", "BHK21"]),
    ("Neuro-2a", "CVCL_0470", ["Neuro 2a", "Neuro2a", "N-2a", "N2a", "Nb2a"]),
    ("Huh-7", "CVCL_0336", ["HuH7", "Huh7", "HUH7.0", "JTC-39"]),
    ("C2C12", "CVCL_0188", ["C2-C12"]),
    ("HEK293-A", "CVCL_6910", ["HEK-293A", "HEK293A", "HEK 293A", "293A", "QBI-293A"]),
    ("HEK293-F", "CVCL_6642", ["HEK-293-F", "HEK 293-F", "HEK-293F", "HEK293F", "293-F"]),
    ("HEK293T/17", "CVCL_1926", ["HEK-293T/17", "HEK 293T/17", "293T/17"]),
    ("SaOS-2", "CVCL_0548", ["SAOS 2", "Saos2"]),
    ("Vero", "CVCL_0059", ["VeroCCL81", "Vero 81", "Vero-81"]),
    ("HT-1080", "CVCL_0317", ["HT 1080", "HT1080"]),
    ("MCF-7", "CVCL_0031", ["MCF 7", "MCF.7", "MCF7"]),
    ("A-549", "CVCL_0023", ["A 549", "A549", "NCI-A549"]),
    ("CV-1", "CVCL_0229", ["CV 1", "CV1"]),
    ("HCT 116", "CVCL_0291", ["HCT-116", "HCT.116", "HCT_116", "HCT116"]),
    # PC12/PC-12\u5728Cellosaurus\u5bf9\u5e94\u4e0d\u540c\u8eab\u4efd\uff1b\u6545\u610f\u4fdd\u7559\u51b2\u7a81\uff0c\u7981\u6b62\u81ea\u52a8\u5408\u5e76\u3002
    ("PC-12", "CVCL_S979", ["PC12", "PC 12"]),
    ("PC12", "CVCL_0481", ["PC-12", "PC 12"]),
]


def text(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def known(value):
    return text(value).strip().casefold() not in UNKNOWN


def protein_type(token):
    body = token.removeprefix("p.").strip("()")
    if "fs" in body:
        return "Frameshift"
    for marker, label in [("delins", "Delins"), ("del", "Deletion"), ("dup", "Duplication"), ("ins", "Insertion")]:
        if marker in body:
            return label
    if re.search(r"(?:Ter|\*)$", body):
        return "Nonsense"
    if body.endswith("="):
        return "Unresolved"  # \u86cb\u767d\u672a\u53d8\u4e0d\u80fd\u81ea\u52a8\u63a8\u51faDNA\u540c\u4e49\u53d8\u5f02\u3002
    m = re.fullmatch("(" + A + r")([1-9]\d*)(" + A + ")", body)
    if m and (m[1] == m[3] or (m[1] == "Met" and m[2] == "1")):
        return "Unresolved"
    return "Missense"


def short_protein(raw, allow_x=False):
    s = raw.strip()
    if s.startswith("p."):
        s = s[2:]
    predicted = s.startswith("(") and s.endswith(")")
    if predicted:
        s = s[1:-1]
    result = None
    m = re.fullmatch(r"([ARNDCEQGHILKMFPSTWYV])([1-9]\d*)([ARNDCEQGHILKMFPSTWYV*X])", s)
    if m and (m[3] != "X" or allow_x):
        result = AA[m[1]] + m[2] + ("Ter" if m[3] in "*X" else AA[m[3]])
    m = re.fullmatch(r"([ARNDCEQGHILKMFPSTWYV])([1-9]\d*)(del|dup)", s)
    if m:
        result = AA[m[1]] + m[2] + m[3]
    m = re.fullmatch(r"[\u0394\u2206]([ARNDCEQGHILKMFPSTWYV])([1-9]\d*)", s)
    if m:
        result = AA[m[1]] + m[2] + "del"
    if re.fullmatch(BODY, s):
        result = s.replace("*", "Ter")
    if result:
        return "p." + ("(" + result + ")" if predicted else result)
    return None


def normalize_mutation(raw, old_type):
    s = text(raw).strip()
    if not known(s):
        return raw, "\u672a\u63d0\u4f9b\u53d8\u5f02\u63cf\u8ff0"
    allow_x = text(old_type).casefold() == "nonsense"
    if P_RE.fullmatch(s):
        return s.replace("*", "Ter"), "\u5df2\u6709HGVS\u5f62\u5f0f\uff1b\u53c2\u8003\u5e8f\u5217\u53ca\u9884\u6d4b/\u5b9e\u6d4b\u72b6\u6001\u672a\u9a8c\u8bc1"
    candidate = short_protein(s, allow_x)
    if candidate:
        return candidate, "\u5b89\u5168\u7b80\u5199\u8f6c\u6362\uff1b\u539f\u7f16\u53f7/\u9884\u6d4b\u62ec\u53f7\u4fdd\u7559\uff0c\u53c2\u8003\u5e8f\u5217\u672a\u9a8c\u8bc1"
    hits = list(dict.fromkeys(m[0].replace("*", "Ter") for m in P_RE.finditer(s)))
    if len(hits) == 1:
        h = hits[0]
        # \u53ea\u5220\u9664\u4e00\u4e2a\u4e25\u683c\u7b49\u4ef7\u7684\u522b\u540d\uff1b\u542b\u9884\u6d4b\u62ec\u53f7\u65f6\u4e0d\u5220\u9664\uff0c\u907f\u514d\u6539\u53d8\u8bc1\u636e\u542b\u4e49\u3002
        if "p.(" not in s:
            for m in re.finditer(r"\([^()]*\)", s):
                inner, outer = m[0][1:-1].strip(), (s[:m.start()] + s[m.end():]).strip()
                if short_protein(inner, allow_x) == h and short_protein(outer, allow_x) == h:
                    return h, "\u5220\u9664\u4e25\u683c\u7b49\u4ef7\u522b\u540d\uff1b\u539f\u6587\u4fdd\u7559\u5728\u539f\u59cb\u6570\u636e/\u5ba1\u8ba1\uff0c\u53c2\u8003\u5e8f\u5217\u672a\u9a8c\u8bc1"
        return raw, "HGVS\u4e0e\u4e0a\u4e0b\u6587\u5e76\u5b58\uff1b\u4fdd\u7559\u539f\u6587\u5f85\u6838"
    if len(hits) > 1:
        return raw, "\u591a\u4e2a\u53d8\u5f02\u6216\u53c2\u7167\u8868\u8fbe\uff1b\u4e0d\u63a8\u65ad\u987a\u53cd\u5173\u7cfb\uff0c\u4fdd\u7559\u539f\u6587"
    return raw, "\u4e0d\u80fd\u5b89\u5168\u8f6c\u6362\uff1a\u4fdd\u7559DNA/RNA\u3001\u6784\u5efa\u4f53\u6216\u4e0d\u5b8c\u6574\u63cf\u8ff0"


def normalize_type(raw, old_type):
    s, t = text(raw), text(old_type).strip().casefold()
    low = s.casefold()
    if s.strip() in {"p.=", "p.(=)"}:
        return "Unresolved", "\u86cb\u767d\u6574\u4f53\u65e0\u53d8\u5316\u4e0d\u7b49\u4e8eDNA\u540c\u4e49\u53d8\u5f02\uff0c\u9700\u6838\u5bf9\u6765\u6e90"
    if t in {"wild-type", "wild type", "wildtype"} or (t == "none" and re.fullmatch(r"wild[- ]type(?: [\w-]+)?;? no mutation", low.strip())):
        return "Wild-type", "\u539f\u7c7b\u578b/\u63cf\u8ff0\u660e\u786e\u91ce\u751f\u578b"
    if t in {"knockdown", "overexpression", "conformational conversion", "non-sequence perturbation"}:
        return "Non-sequence perturbation", "\u975e\u5e8f\u5217\u6270\u52a8\u5355\u5217"
    # \u53ea\u5bf9\u5b8c\u6574\u8868\u8fbe\u5f0f\u6216\u4e25\u683c\u7b49\u4ef7\u522b\u540d\u5224\u65ad\uff1b\u53e5\u4e2d\u7684p.\u53ef\u80fd\u662f\u5bf9\u7167/\u5176\u4ed6\u6784\u5efa\u4f53\u3002
    cleaned, _ = normalize_mutation(raw, old_type)
    simple = short_protein(text(cleaned), t == "nonsense")
    hits = [simple] if simple else []
    categories = {protein_type(h) for h in hits}
    if len(categories) == 1 and "Unresolved" not in categories:
        return next(iter(categories)), "\u660e\u786e\u86cb\u767d\u540e\u679c\u4f18\u5148\uff1bfs\u5148\u4e8e\u7ec8\u6b62\uff0cdelins\u5148\u4e8edel/ins"
    if categories == {"Unresolved"}:
        return ("Synonymous", "\u4fdd\u7559\u660e\u786e\u539f\u7c7b\u578b\uff1bp.=\u672c\u8eab\u4e0d\u8bc1\u660eDNA\u540c\u4e49\u53d8\u5f02") if t == "synonymous" else ("Unresolved", "\u86cb\u767d\u672a\u53d8\u6216\u8d77\u59cb\u4f4d\u70b9\u540e\u679c\u9700\u6838\u5b9e")
    if re.fullmatch(r"frame[- ]?shift(?: mutation)?", low.strip()):
        return "Frameshift", "\u72ec\u7acb\u660e\u786e\u7684\u79fb\u7801\u4e8b\u4ef6\u540d\u79f0"
    contextual = bool(P_RE.search(s) or re.search(r"frame[- ]?shift", low)) and not simple
    if t in TYPE_ALIASES:
        return TYPE_ALIASES[t], ("\u4e0a\u4e0b\u6587/\u6bd4\u8f83/\u5426\u5b9a\u6216\u4e0d\u786e\u5b9a\u63cf\u8ff0\uff0c\u4e0d\u636e\u6b64\u8986\u76d6\u539f\u7c7b\u578b\uff1b\u5f85\u6838" if contextual else "\u5df2\u6709\u7c7b\u578b\u540d\u79f0\u5f52\u5e76\uff1b\u672a\u9a8c\u8bc1\u5e8f\u5217")
    if not known(t) or t in {"mutation", "point mutation", "polymorphism", "natural variant", "none", "null", "hypomorphic", "loss-of-function"}:
        return "Unresolved", "\u4fe1\u606f\u4e0d\u8db3\uff0c\u4e0d\u5f3a\u5236\u63a8\u65ad\u53d8\u5f02\u7c7b\u522b"
    return "Other", "\u5176\u4ed6\u7ed3\u6784\u53d8\u5f02/\u622a\u77ed/\u6784\u5efa\u4f53\uff0c\u539f\u7c7b\u578b\u4fdd\u7559"


def cell_key(s):
    return re.sub(r"[\s\-\u2010\u2011\u2212]", "", s).casefold()


def build_cell_index():
    index = defaultdict(dict)
    for name, accession, aliases in CELL_LINES:
        for alias in [name, *aliases]:
            index[cell_key(alias)][accession] = name
    return index


CELL_INDEX = build_cell_index()


def normalize_cell(value):
    s = text(value)
    if not known(s):
        return value, [], ["\u7ec6\u80de\u540d\u79f0\u7f3a\u5931/\u672a\u786e\u5b9a"]
    # \u5b8c\u6574\u540c\u4e49\u8bcd\u4f18\u5148\uff0c\u907f\u514d\u628a HEK;293 \u8bef\u5f53\u4e24\u4e2a\u7ec6\u80de\u3002
    tokens = [s] if cell_key(s.strip()) in CELL_INDEX else s.split(";")
    outputs, accessions, flags = [], [], []
    for token in tokens:
        token = token.strip()
        key = cell_key(token)
        matches = CELL_INDEX.get(key)
        if matches is None:
            key = cell_key(re.sub(r"\s+cells?$", "", token, flags=re.I))
            matches = CELL_INDEX.get(key, {})
        if len(matches) == 1:
            accession, name = next(iter(matches.items()))
            outputs.append(name)
            accessions.append(accession)
        else:
            outputs.append(token)
            flags.append("\u7ec6\u80de\u540d\u79f0\u6b67\u4e49\uff0c\u4fdd\u7559\u539f\u6587" if matches else "\u7ec6\u80de\u540d\u79f0\u672a\u5339\u914d\uff0c\u4fdd\u7559\u539f\u6587/\u6761\u4ef6")
    return "; ".join(outputs), list(dict.fromkeys(accessions)), list(dict.fromkeys(flags))


def header_key(s):
    return re.sub(r"[\s_()\uff08\uff09\-]", "", text(s)).casefold()


HEADER_LOOKUP = {header_key(h): h for h in BASE_COLUMNS + [DETAIL]}
HEADER_LOOKUP.update({header_key(k): v for k, v in {
    "gene_name": "Gene Name", "mutation_information": "Mutation Information",
    "normal_localization": "Mislocalization Phenotype\uff08before\uff09",
    "mutated_localization": "Mislocalization Phenotype\uff08after\uff09",
    "literature_title": "Title", "source_pmid": "PMID", "record_id": "id",
}.items()})


def canonical_header(value):
    s = text(value).strip()
    k = header_key(s)
    if k.startswith("mislocalizationphenotype"):
        if "before" in k:
            return "Mislocalization Phenotype\uff08before\uff09"
        if "after" in k:
            return "Mislocalization Phenotype\uff08after\uff09"
    return HEADER_LOOKUP.get(k, s)


def sheet_header(ws):
    return [canonical_header(v) for v in next(ws.iter_rows(values_only=True), ())]


def select_sheet(wb, requested):
    if requested:
        if requested not in wb.sheetnames:
            raise ValueError(f"\u6ca1\u6709\u5de5\u4f5c\u8868 {requested!r}\uff1b\u53ef\u7528\uff1a{', '.join(wb.sheetnames)}")
        candidates = [requested]
    else:
        candidates = [s.title for s in wb if REQUIRED <= set(sheet_header(s))]
        for preferred in ["\u5b8c\u6574\u8bb0\u5f55_\u65e7\u8868", "5776\u6761\u89c4\u8303\u8868", "\u89c4\u8303\u5316\u6570\u636e"]:
            if preferred in candidates:
                return wb[preferred]
        if len(candidates) != 1:
            raise ValueError(f"\u65e0\u6cd5\u552f\u4e00\u9009\u62e9\u6570\u636e\u8868\uff0c\u8bf7\u7528 --sheet \u6307\u5b9a\u3002\u5019\u9009\uff1a{candidates}")
    if not candidates or not REQUIRED <= set(sheet_header(wb[candidates[0]])):
        raise ValueError("\u6570\u636e\u8868\u81f3\u5c11\u9700\u8981 PMID\u3001Mutation Information\u3001Mislocalization Phenotype \u4e09\u5217")
    return wb[candidates[0]]


def iter_table(ws):
    """\u6821\u9a8c\u8868\u5934\u540e\u8fd4\u56de\u884c\u8fed\u4ee3\u5668\uff0c\u8f85\u52a9\u8868\u4e0d\u6574\u4f53\u7269\u5316\u5230\u5185\u5b58\u3002"""
    iterator = ws.iter_rows(values_only=True)
    headers = [canonical_header(v) for v in next(iterator, ())]
    while headers and not headers[-1]:
        headers.pop()
    if not headers or any(not h for h in headers) or len(set(headers)) != len(headers):
        raise ValueError(f"\u5de5\u4f5c\u8868 {ws.title} \u7684\u8868\u5934\u4e3a\u7a7a\u6216\u89c4\u8303\u540e\u91cd\u540d\uff0c\u8bf7\u5148\u4fee\u6b63")
    def rows():
        for excel_row, cells in enumerate(iterator, 2):
            if not any(v is not None for v in cells):
                continue
            if any(v is not None for v in cells[len(headers):]):
                raise ValueError(f"{ws.title} \u7b2c{excel_row}\u884c\u5b58\u5728\u65e0\u8868\u5934\u6570\u636e")
            yield excel_row, dict(zip(headers, list(cells[:len(headers)]) + [None] * max(0, len(headers) - len(cells))))
    return headers, rows()


def read_table(ws):
    headers, iterator = iter_table(ws)
    rows = list(iterator)
    if not rows:
        raise ValueError(f"\u5de5\u4f5c\u8868 {ws.title} \u6ca1\u6709\u8bb0\u5f55")
    return headers, rows


def auxiliary_index(wb, selected, wanted_keys=None):
    """\u4ec5\u6309\u53cc\u952e\u8fde\u63a5\u540c\u5de5\u4f5c\u7c3f\u8f85\u52a9\u8868\uff1b\u53cc\u952e\u51b2\u7a81\u4e0d\u91c7\u7528\uff0c\u9632\u6b62\u9519\u7ed1\u3002"""
    index, conflicts = {}, set()
    for name in ["\u590d\u6838\u7d22\u5f15", "\u89c4\u8303\u5ba1\u8ba1", "\u5904\u7406\u72b6\u6001"]:
        if name not in wb.sheetnames or name == selected:
            continue
        ws = wb[name]
        if not {"id", "PMID"} <= set(sheet_header(ws)):
            continue
        _, rows = iter_table(ws)
        seen = set()
        for _, row in rows:
            key = (text(row.get("id")).strip(), text(row.get("PMID")).strip())
            if not all(key) or (wanted_keys is not None and key not in wanted_keys):
                continue
            if name == "\u5904\u7406\u72b6\u6001":
                row = {"__prior_flags": row.get("\u5f85\u6838\u4e8b\u9879"), "__prior_detail_basis": row.get("\u5b9a\u4f4d\u7ec6\u8282\u4f9d\u636e")}
            if key in seen:
                conflicts.add(key)
            seen.add(key)
            target = index.setdefault(key, {})
            for field, value in row.items():
                if field not in {"id", "PMID"} and value is not None:
                    if field in target and target[field] != value:
                        conflicts.add(key)
                    target[field] = value
    for key in conflicts:
        index.pop(key, None)
    return index, conflicts


def field(row, *names):
    for name in names:
        if known(row.get(name)):
            return text(row[name])
    return ""


def localization_detail(row, aux):
    facts = dict(aux)
    facts.update({k: v for k, v in row.items() if known(v)})
    before = field(row, "Mislocalization Phenotype\uff08before\uff09")
    after = field(row, "Mislocalization Phenotype\uff08after\uff09")
    same = bool(before and after and before.strip().casefold() == after.strip().casefold())
    different = bool(before and after and not same)
    outcome = field(row, "Outcome_Provenance").strip().upper()
    event = field(facts, "localization_event_type", "Localization Event Type").strip().casefold()
    extent = field(facts, "localization_extent").strip().casefold()
    review_extent = field(row, "Redistribution_Extent").strip().casefold()
    flags = []
    names = {"retention": "\u6ede\u7559", "aggregation": "\u805a\u96c6", "fragmentation": "\u788e\u7247\u5316",
             "redistribution": "\u518d\u5206\u5e03\uff08\u5177\u4f53\u5f62\u5f0f\u5f85\u6838\uff09", "translocation": "\u8de8\u533a\u5ba4\u8f6c\u4f4d", "no_change": "\u5b9e\u9a8c\u65e0\u53d8\u5316"}
    changing = set(names) - {"no_change"}
    no_change = outcome == "EXPERIMENTAL_NO_CHANGE" or event == "no_change"
    extents = {v for v in (extent, review_extent) if known(v) and v != "not_applicable"}
    change_extents = extents - {"no_change", "none", "insufficient"}
    if same and event == "translocation":
        flags.append("\u540c\u533a\u5ba4\u7f16\u7801\u4e0e\u8de8\u533a\u5ba4\u8f6c\u4f4d\u4e8b\u4ef6\u51b2\u7a81")
    if no_change and (different or event in changing or change_extents):
        flags.append("\u65e0\u53d8\u5316\u6807\u7b7e\u4e0e\u533a\u5ba4/\u4e8b\u4ef6/\u53d8\u5316\u7a0b\u5ea6\u51b2\u7a81")
    if outcome == "EXPERIMENTAL_MISLOCALIZATION" and (event == "no_change" or "no_change" in extents):
        flags.append("\u5b9e\u9a8c\u9519\u5b9a\u4f4d\u4e0e\u65e0\u53d8\u5316\u6807\u7b7e\u51b2\u7a81")
    if "no_change" in extents and (different or event in changing or change_extents):
        flags.append("\u65e0\u53d8\u5316\u7a0b\u5ea6\u4e0e\u533a\u5ba4/\u4e8b\u4ef6\u51b2\u7a81")
    existing = field(row, DETAIL)
    if existing:
        if no_change and (existing in {"\u6ede\u7559", "\u805a\u96c6", "\u788e\u7247\u5316", "\u8de8\u533a\u5ba4\u8f6c\u4f4d", "\u5b9a\u4f4d\u6bd4\u4f8b/\u5b9a\u91cf\u518d\u5206\u5e03"} or existing.startswith("\u518d\u5206\u5e03")):
            flags.append("\u5df2\u6709\u5b9a\u4f4d\u7ec6\u8282\u4e0e\u65e0\u53d8\u5316\u6807\u7b7e\u51b2\u7a81")
        if outcome == "EXPERIMENTAL_MISLOCALIZATION" and existing in {"\u5b9e\u9a8c\u65e0\u53d8\u5316", "\u6e90\u8bb0\u5f55\u6807\u4e3a\u5b9e\u9a8c\u65e0\u53d8\u5316"}:
            flags.append("\u5df2\u6709\u65e0\u53d8\u5316\u7ec6\u8282\u4e0e\u5b9e\u9a8c\u9519\u5b9a\u4f4d\u51b2\u7a81")
        if "\u5f85\u660e\u786e" in existing:
            flags.append("\u673a\u5236\u4fe1\u606f\u4e0d\u8db3")
        if any(word in existing for word in ["\u5f85\u6838", "\u5f85\u590d\u6838", "\u63d0\u793a"]):
            flags.append("\u5b9a\u4f4d\u7ec6\u8282\u9700\u539f\u6587\u786e\u8ba4")
        if "\u51b2\u7a81" in existing:
            flags.append("\u5b9a\u4f4d\u673a\u5236\u6807\u7b7e\u51b2\u7a81\u5f85\u590d\u6838")
        return row[DETAIL], field(aux, "__prior_detail_basis") or "\u4fdd\u7559\u8f93\u5165\u5df2\u6709\u5b9a\u4f4d\u7ec6\u8282\u5e76\u91cd\u65b0\u68c0\u67e5", flags
    if flags:
        return "\u5b9a\u4f4d\u533a\u5ba4/\u4e8b\u4ef6/\u7ed3\u679c\u6807\u7b7e\u51b2\u7a81\uff0c\u5f85\u590d\u6838", f"event={event}; outcome={outcome}; extent={extent}/{review_extent}", flags
    if event in names:
        detail = names[event]
        if event == "redistribution" and extent == "quantitative":
            detail = "\u5b9a\u4f4d\u6bd4\u4f8b/\u5b9a\u91cf\u518d\u5206\u5e03"
        if "\u5f85\u6838" in detail:
            flags.append("\u5b9a\u4f4d\u7ec6\u8282\u9700\u539f\u6587\u786e\u8ba4")
        return detail, f"\u8868\u5185\u7ed3\u6784\u5316\u4e8b\u4ef6\u7c7b\u578b={event}\uff1b\u4e0d\u4ee3\u8868\u4eba\u5de5\u6838\u9a8c", flags
    if outcome == "EXPERIMENTAL_NO_CHANGE":
        return "\u6e90\u8bb0\u5f55\u6807\u4e3a\u5b9e\u9a8c\u65e0\u53d8\u5316", "\u8868\u5185Outcome_Provenance", flags
    # \u539f\u59cb\u5b9a\u4f4d\u63cf\u8ff0\u4ec5\u751f\u6210\u63d0\u793a\u3002\u8bc1\u636e\u539f\u6587\u53ef\u80fd\u6d89\u53ca\u5176\u4ed6\u86cb\u767d\uff0c\u6545\u4e0d\u505a\u5173\u952e\u8bcd\u81ea\u52a8\u5b9a\u8bba\u3002
    raw = field(facts, "mutated_localization_raw", "\u539fafter\u63cf\u8ff0", "\u539f\u59cb\u7a81\u53d8\u5b9a\u4f4d")
    candidates = []
    for pattern, label in [(r"\bret(?:ain|ention)\w*|\u6ede\u7559", "\u6ede\u7559"),
                           (r"\baggregat\w*|\u805a\u96c6", "\u805a\u96c6"),
                           (r"\bfragment\w*|\u788e\u7247\u5316", "\u788e\u7247\u5316"),
                           (r"\bco[- ]?locali[sz]ation\b|\u5171\u5b9a\u4f4d", "\u5171\u5b9a\u4f4d\u53d8\u5316")]:
        if raw and re.search(pattern, raw, re.I):
            if not re.search(r"\b(?:no|not|without|absence)\b|\u6ca1\u6709|\u672a\u89c1|\u65e0", raw, re.I):
                candidates.append(label)
    if candidates:
        return "\u539f\u59cb\u5b9a\u4f4d\u63cf\u8ff0\u63d0\u793a" + "/".join(candidates) + "\uff08\u5f85\u590d\u6838\uff09", "\u540c\u5de5\u4f5c\u7c3f\u539f\u59cb\u5b9a\u4f4d\u63cf\u8ff0\uff0c\u975e\u673a\u5236\u786e\u8ba4", ["\u5b9a\u4f4d\u7ec6\u8282\u9700\u539f\u6587\u786e\u8ba4"]
    return ("\u540c\u533a\u5ba4\u53d8\u5316\uff0c\u673a\u5236\u5f85\u660e\u786e" if same else "\u673a\u5236\u5f85\u660e\u786e"), "XLSX\u7f3a\u5c11\u8db3\u591f\u673a\u5236\u4fe1\u606f\uff0c\u4e0d\u4ece\u533a\u5ba4\u7f16\u7801\u731c\u6d4b", ["\u673a\u5236\u4fe1\u606f\u4e0d\u8db3"]


class TitleResolver:
    def __init__(self, mode="online", cache_path=None, timeout=20, retries=1, refresh=False, max_age_days=None):
        self.mode, self.cache_path, self.timeout, self.retries = mode, cache_path, timeout, retries
        self.refresh, self.max_age_days = refresh, max_age_days
        self.refreshed = set()
        self.entries, self.errors = {}, []
        if cache_path and cache_path.exists():
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
            if not isinstance(cache, dict) or cache.get("schema") != "xlsx-normalizer-title-cache-v1":
                raise ValueError("\u6807\u9898\u7f13\u5b58\u683c\u5f0f\u4e0d\u5339\u914d\uff1b\u8bf7\u4f7f\u7528\u672c\u7a0b\u5e8f\u4ea7\u751f\u7684\u7f13\u5b58\u6216\u66f4\u6362 --cache")
            self.entries = cache.get("titles", {})
            if not isinstance(self.entries, dict):
                raise ValueError("\u6807\u9898\u7f13\u5b58 titles \u5fc5\u987b\u662f\u5bf9\u8c61")

    def get_json(self, base, params):
        url = base + "?" + urlencode(params)
        for attempt in range(self.retries + 1):
            try:
                request = Request(url, headers={"User-Agent": "MutProML-XLSX-Normalizer/" + VERSION,
                                                "Accept": "application/json"})
                with urlopen(request, timeout=self.timeout) as response:
                    data = json.load(response)
                if not isinstance(data, dict):
                    raise ValueError("API\u54cd\u5e94\u4e0d\u662fJSON\u5bf9\u8c61")
                return data
            except (OSError, ValueError) as exc:
                if attempt == self.retries:
                    raise RuntimeError(f"\u5143\u6570\u636e\u8bf7\u6c42\u5931\u8d25 {base}: {exc}") from exc
                time.sleep(min(2 ** attempt, 4))

    def remember(self, pmid, result, method, url, binding_url=""):
        title = result.get("title")
        if not isinstance(title, str) or not title.strip():
            return
        self.entries[pmid] = {"pmid": pmid, "title": html.unescape(re.sub(r"<[^>]+>", "", title)).strip(),
                              "method": method, "url": url, "binding_url": binding_url,
                              "retrieved_at": datetime.now(timezone.utc).isoformat()}
        self.refreshed.add(pmid)

    def checkpoint(self):
        if self.cache_path:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            atomic_json(self.cache_path, {"schema": "xlsx-normalizer-title-cache-v1", "titles": self.entries})

    def usable(self, pmid):
        e = self.entries.get(pmid)
        valid = isinstance(e, dict) and e.get("pmid") == pmid and isinstance(e.get("title"), str) and bool(e["title"].strip()) and bool(e.get("url"))
        if not valid or (self.refresh and self.mode == "online" and pmid not in self.refreshed):
            return False
        if self.max_age_days is not None:
            try:
                retrieved = datetime.fromisoformat(e.get("retrieved_at", ""))
                age = (datetime.now(timezone.utc) - retrieved).total_seconds() / 86400
                if age < 0 or age > self.max_age_days:
                    return False
            except (ValueError, TypeError):
                return False
        return True

    def resolve(self, pmids):
        valid = sorted({p for p in pmids if re.fullmatch(r"[1-9]\d*", p)}, key=int)
        if self.mode != "online":
            return
        missing = [p for p in valid if not self.usable(p)]
        failed_streak = 0
        for start in range(0, len(missing), 80):
            batch = missing[start:start + 80]
            query = "SRC:MED AND (" + " OR ".join("EXT_ID:" + p for p in batch) + ")"
            try:
                data = self.get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                                     {"query": query, "format": "json", "pageSize": 1000, "resultType": "lite"})
                for result in data.get("resultList", {}).get("result", []):
                    pmid = text(result.get("pmid") or result.get("id"))
                    if result.get("source") == "MED" and pmid in batch and text(result.get("id")) == pmid:
                        self.remember(pmid, result, "Europe PMC MED PMID\u7cbe\u786e\u6838\u5b9a", "https://europepmc.org/article/MED/" + pmid)
                failed_streak = 0
            except RuntimeError as exc:
                self.errors.append(str(exc))
                failed_streak += 1
            self.checkpoint()  # \u6bcf\u6279\u539f\u5b50\u843d\u76d8\uff0c\u540e\u7eed\u4e2d\u65ad\u65f6\u5df2\u6709\u6279\u6b21\u53ef\u590d\u7528\u3002
            print(f"\u6807\u9898\u6838\u5b9a\uff1a{sum(self.usable(p) for p in valid)}/{len(valid)}", flush=True)
            if failed_streak >= 3:
                self.errors.append("\u8fde\u7eed\u4e09\u6279\u7f51\u7edc\u5931\u8d25\uff0c\u505c\u6b62\u672c\u6b21\u540e\u7eed\u8bf7\u6c42\uff1b\u672a\u6838\u5b9a\u6807\u9898\u4fdd\u7559\u539f\u503c")
                break
            time.sleep(0.35)
        # \u7528PMC\u5b98\u65b9\u6620\u5c04\u6838\u5b9a\u9884\u5370\u672c/\u672a\u547d\u4e2dMED\u7684PMID\uff1b\u4e0d\u4f9d\u8d56XLSX\u91cc\u7684\u672a\u6838\u5b9aDOI\u3002
        rest = [p for p in valid if not self.usable(p)]
        if failed_streak < 3:
            for start in range(0, len(rest), 200):
                batch = rest[start:start + 200]
                try:
                    bindings = self.get_json("https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/",
                                             {"ids": ",".join(batch), "format": "json", "tool": "MutProML_XLSX_Normalizer"})
                    for binding in bindings.get("records", []):
                        pmid, doi = text(binding.get("pmid")), binding.get("doi")
                        if pmid not in batch or not isinstance(doi, str) or not doi:
                            continue
                        try:
                            data = self.get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                                                 {"query": 'DOI:"' + doi.replace('"', '') + '"', "format": "json", "resultType": "lite"})
                            hits = [r for r in data.get("resultList", {}).get("result", []) if text(r.get("doi")).casefold() == doi.casefold()]
                            if len(hits) == 1:
                                hit = hits[0]
                                self.remember(pmid, hit, "PMC PMID\u2192DOI\u2192Europe PMC\u7cbe\u786e\u6838\u5b9a\uff08\u4fdd\u7559\u539f\u7248\u672c\uff09",
                                              "https://europepmc.org/article/" + text(hit.get("source")) + "/" + text(hit.get("id")),
                                              "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?ids=" + pmid + "&format=json")
                        except RuntimeError as exc:
                            self.errors.append(str(exc))
                        self.checkpoint()  # DOI\u56de\u9000\u6bcf\u6761\u6838\u5b9a\u540e\u4fdd\u5b58\u3002
                        time.sleep(0.35)
                except RuntimeError as exc:
                    self.errors.append(str(exc))
        self.checkpoint()


def atomic_json(path, value):
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def parse_key_values(value):
    """\u53ea\u7528\u4e8e\u68c0\u67e5\uff1b\u4e0d\u6539\u5199\u539f\u5b57\u6bb5\uff0c\u91cd\u590d\u952e\u4e5f\u4e0d\u9759\u9ed8\u8986\u76d6\u3002"""
    result = defaultdict(set)
    for part in re.split(r"[;\uff1b]", text(value)):
        if "=" in part:
            key, val = part.split("=", 1)
            result[key.strip().casefold()].add(val.strip().upper())
    return result


def review_flags(row):
    flags = []
    tag = parse_key_values(row.get("Tag_Strategy"))
    if any(len(values) > 1 for values in tag.values()):
        flags.append("\u6807\u7b7e\u7b56\u7565\u542b\u51b2\u7a81\u7684\u91cd\u590d\u952e\uff0c\u9700\u67e5\u539f\u6587")
    if tag.get("kind", set()) & {"FLUORESCENT_PROTEIN", "EPITOPE_AFFINITY"} and "NONE" in tag.get("position", set()):
        flags.append("\u6709\u6807\u7b7e\u4f46\u4f4d\u7f6eNONE\uff0c\u9700\u67e5\u539f\u6587\uff1b\u672a\u81ea\u52a8\u6539\u5199")
    provenance, evidence = text(row.get("Outcome_Provenance")).strip().upper(), text(row.get("Evidence_Type")).strip().upper()
    if provenance.startswith("EXPERIMENTAL_") and evidence in {"FUNCTION_ONLY", "SECONDARY_ONLY", "UNDETERMINED"}:
        flags.append("\u5b9e\u9a8c\u7ed3\u679c\u6807\u7b7e\u4e0e\u8bc1\u636e\u7c7b\u578b\u9700\u8054\u5408\u590d\u6838")
    if provenance == "EXPERIMENTAL_NO_CHANGE" and text(row.get("No_Change_Basis")).strip().upper() in {"", "INSUFFICIENT", "NOT_APPLICABLE", "UNKNOWN", "UNDETERMINED", "NOT_REPORTED"}:
        flags.append("\u5b9e\u9a8c\u65e0\u53d8\u5316\u7684\u5bf9\u7167\u4f9d\u636e\u5f85\u6838")
    if provenance == "EXPERIMENTAL_NO_CHANGE" and text(row.get("Redistribution_Extent")).strip().upper() not in {"", "NO_CHANGE"}:
        flags.append("\u65e0\u53d8\u5316\u6807\u7b7e\u4e0e\u53d8\u5316\u7a0b\u5ea6\u4e0d\u4e00\u81f4")
    if any(isinstance(v, str) and v.startswith("=") for v in row.values()):
        flags.append("\u8f93\u5165\u516c\u5f0f\u4f5c\u4e3a\u539f\u59cb\u6587\u672c\u4fdd\u7559\uff0c\u672a\u6267\u884c/\u8ba1\u7b97")
    return flags


def normalize_rows(headers, source_rows, aux_index, conflicts, resolver):
    out_headers = BASE_COLUMNS[:6] + [DETAIL] + BASE_COLUMNS[6:]
    out_headers += [h for h in headers if h not in out_headers]
    normalized, changes, statuses = [], [], []
    mappings, counts = Counter(), Counter()
    seen_ids = Counter(text(r.get("id")).strip() for _, r in source_rows if r.get("id") is not None)
    seen_keys = Counter((text(r.get("id")).strip(), text(r.get("PMID")).strip()) for _, r in source_rows if known(r.get("id")))
    occupied_ids = set(seen_ids)
    for serial, (excel_row, original) in enumerate(source_rows, 1):
        row = {h: original.get(h) for h in out_headers}
        had_id = known(row.get("id"))
        if not had_id:
            generated = f"xlsx-row-{excel_row}"
            while generated in occupied_ids:
                generated += "_"
            row["id"] = generated
            occupied_ids.add(generated)
        row["PMID"] = text(row.get("PMID")).strip()
        rid, pmid = text(row["id"]), row["PMID"]
        key = (rid.strip(), pmid)
        flags = review_flags(original)
        if not had_id:
            flags.append("\u7f3a\u5c11id\uff0c\u6309Excel\u884c\u751f\u6210\u552f\u4e00\u6807\u8bc6\uff1b\u672a\u8fde\u63a5\u8f85\u52a9\u8868")
        if seen_ids.get(rid.strip(), 0) > 1:
            flags.append("\u8f93\u5165id\u91cd\u590d\uff1b\u4f7f\u7528\u6765\u6e90Excel\u884c\u8ffd\u6eaf\uff0c\u4e0d\u53bb\u91cd")
        if seen_keys[key] > 1:
            flags.append("\u4e3b\u8868\u53cc\u952e\u4e0d\u552f\u4e00\uff0c\u672a\u8fde\u63a5\u8f85\u52a9\u8bc1\u636e")
            counts["main_duplicate_key_rows"] += 1
        if key in conflicts:
            flags.append("\u8f85\u52a9\u8868\u53cc\u952e\u4e0d\u552f\u4e00/\u51b2\u7a81\uff0c\u672a\u91c7\u7528\u5176\u5185\u5bb9")
            counts["aux_conflict_rows"] += 1
        aux = aux_index.get(key, {}) if had_id and seen_keys[key] == 1 and key not in conflicts else {}
        # \u4ec5\u7ee7\u627f\u53ef\u9760\u53cc\u952e\u8fde\u63a5\u7684\u590d\u6838\u8bb0\u5f55\uff1b\u6807\u9898\u72b6\u6001\u6309\u672c\u6b21\u6838\u5b9a\u91cd\u65b0\u8ba1\u7b97\u3002
        prior = text(aux.get("__prior_flags"))
        flags.extend(f for f in prior.split("; ") if f and f != "\u6807\u9898\u5f85\u5916\u90e8\u6838\u5b9a")
        if not re.fullmatch(r"[1-9]\d*", pmid):
            flags.append("PMID\u683c\u5f0f\u65e0\u6548\uff0c\u672a\u8054\u7f51\u67e5\u8be2")
        variant, hgvs_status = normalize_mutation(row.get("Mutation Information"), row.get("Mutation Type"))
        mtype, type_reason = normalize_type(row.get("Mutation Information"), row.get("Mutation Type"))
        if "\u5f85\u6838" in type_reason:
            flags.append("\u53d8\u5f02\u63cf\u8ff0\u542b\u4e0a\u4e0b\u6587/\u5426\u5b9a/\u6bd4\u8f83\uff0c\u7c7b\u578b\u5f85\u6838")
        row["Mutation Information"], row["Mutation Type"] = variant, mtype
        mappings[(text(original.get("Mutation Type")), mtype, type_reason)] += 1
        row["Cell Line"], rrids, cell_flags = normalize_cell(row.get("Cell Line"))
        flags.extend(cell_flags)
        row[DETAIL], detail_basis, detail_flags = localization_detail(row, aux)
        flags.extend(detail_flags)
        for name in BLANK_COLUMNS:
            row[name] = None
        if resolver.mode != "keep" and resolver.usable(pmid):
            row["Title"] = resolver.entries[pmid]["title"]
            title_status = resolver.entries[pmid].get("method", "\u8bfb\u53d6\u5df2\u6838\u5b9a\u7f13\u5b58")
        else:
            title_status = "\u4fdd\u7559\u539f\u6807\u9898\uff0c\u672a\u5b8c\u6210\u5916\u90e8\u6838\u5b9a"
            flags.append("\u6807\u9898\u5f85\u5916\u90e8\u6838\u5b9a")
            counts["title_unverified_rows"] += 1
        if mtype in {"Other", "Complex", "Unresolved", "Wild-type", "Non-sequence perturbation"}:
            flags.append("\u53d8\u5f02/\u6784\u5efa\u4f53\u7c7b\u522b\u9700\u7ed3\u5408\u7eb3\u5165\u8303\u56f4\u590d\u6838")
        if "p." in text(variant):
            flags.append("HGVS\u53c2\u8003\u5e8f\u5217/\u7f16\u53f7/\u9884\u6d4b\u4e0e\u5b9e\u6d4b\u72b6\u6001\u672a\u9a8c\u8bc1")
        row_changed = False
        for name in out_headers:
            old, new = original.get(name), row[name]
            if old != new and not (old in (None, "") and new in (None, "")):
                why = type_reason if name == "Mutation Type" else hgvs_status if name == "Mutation Information" else detail_basis if name == DETAIL else title_status if name == "Title" else "\u6309\u7528\u6237\u8981\u6c42\u7559\u7a7a\uff1b\u539f\u503c\u4fdd\u7559" if name in BLANK_COLUMNS else "Cellosaurus\u4fdd\u5b88\u540c\u4e49\u8bcd\u6620\u5c04/\u5206\u9694\u7b26\u6574\u7406" if name == "Cell Line" else "\u6807\u8bc6\u7b26/\u5217\u7ed3\u6784\u89c4\u8303"
                changes.append([excel_row, rid, pmid, name, old, new, why])
                counts["changed_" + name] += 1
                row_changed = True
        assert all(row.get(h) == original.get(h) for h in ["Gene Name", "Disease Information", *BASE_COLUMNS[16:]])
        assert all(row[h] is None for h in BLANK_COLUMNS)
        statuses.append([excel_row, rid, pmid, hgvs_status, detail_basis, "; ".join(rrids), title_status, "; ".join(dict.fromkeys(flags))])
        counts["changed_rows"] += int(row_changed)
        counts["pending_rows"] += bool(flags)
        counts["hgvs_unconverted_rows"] += int(not (P_RE.fullmatch(text(variant)) or text(variant) in {"p.=", "p.(=)"}))
        counts["cell_unmatched_rows"] += int(bool(cell_flags))
        counts["mechanism_pending_rows"] += int(any("\u673a\u5236" in f or "\u5b9a\u4f4d\u7ec6\u8282" in f or "\u51b2\u7a81" in f for f in flags))
        normalized.append(row)
    counts["records"] = len(normalized)
    counts["columns"] = len(out_headers)
    counts["unique_pmids"] = len({r["PMID"] for r in normalized})
    counts["main_duplicate_keys"] = sum(n > 1 for n in seen_keys.values())
    return out_headers, normalized, changes, statuses, mappings, counts


def excel_value(value):
    if isinstance(value, (dict, list, tuple)):
        value = json.dumps(value, ensure_ascii=False)
    if isinstance(value, str):
        value = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\ufffe\uffff]", lambda m: "\\u" + format(ord(m[0]), "04x"), value)
        if len(value) > 32767:
            raise ValueError("\u5355\u5143\u683c\u8d85\u8fc7Excel\u768432767\u5b57\u7b26\u9650\u5236\uff0c\u5df2\u505c\u6b62\u800c\u975e\u622a\u65ad\u6570\u636e")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("\u8f93\u5165\u542bNaN/Infinity\uff0c\u4e0d\u80fd\u65e0\u635f\u5199\u5165Excel")
    return value


def update_digest(digest, row):
    digest.update(json.dumps([None if v in (None, "") else excel_value(v) for v in row], ensure_ascii=False, default=str, separators=(",", ":")).encode("utf-8"))


def data_digest(rows):
    digest = hashlib.sha256()
    for row in rows:
        update_digest(digest, row)
    return digest.hexdigest()


def write_sheet(wb, name, headers, rows, widths=None, *, plain=False):
    from openpyxl.cell import WriteOnlyCell
    ws = wb.create_sheet(name)
    if plain:
        from openpyxl.utils import get_column_letter
        for index, header in enumerate(headers, 1):
            ws.column_dimensions[get_column_letter(index)].width = 22 if index == 1 else min(36, max(18, len(header) + 2))
        digest = hashlib.sha256()
        for values in chain([headers], rows):
            if len(values) != len(headers):
                raise ValueError(f'{name}: wrong column count')
            update_digest(digest, values)
            cells = []
            for value in values:
                # Explicit empty cells keep all 26 positions in streaming XLSX
                # readback without requiring a style to materialize the cell.
                cell = WriteOnlyCell(ws, value='' if value is None else excel_value(value))
                if isinstance(cell.value, str):
                    cell.data_type = 's'
                cells.append(cell)
            ws.append(cells)
        return digest.hexdigest()
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "C2" if len(headers) > 8 else "A2"
    widths = widths or [min(55, max(20, len(h) + 3)) for h in headers]
    for i, width in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = width
    base_font, header_font = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True, color="FFFFFF")
    fill = PatternFill("solid", fgColor="203E5B")
    align = Alignment(vertical="top", wrap_text=True)
    digest = hashlib.sha256()
    for row_number, values in enumerate(chain([headers], rows), 1):
        update_digest(digest, values)
        if len(values) != len(headers):
            raise ValueError(f"{name} \u884c{row_number}\u5217\u6570\u9519\u8bef")
        lengths = [sum(max(1, math.ceil(len(t) / max(8, widths[k] - 5))) for t in text(v).split("\n")) for k, v in enumerate(values)]
        ws.row_dimensions[row_number].height = 48 if row_number == 1 else min(409, max(32, max(lengths, default=1) * 15 + 10))
        cells = []
        for v in values:
            c = WriteOnlyCell(ws, value=excel_value(v))
            if isinstance(c.value, str):
                c.data_type = "s"  # \u5305\u62ec=\u5f00\u5934\u7684\u539f\u6587\uff0c\u7981\u6b62\u88abExcel\u5f53\u4f5c\u516c\u5f0f\u3002
            c.font, c.alignment = (header_font if row_number == 1 else base_font), align
            if row_number == 1:
                c.fill = fill
            cells.append(c)
        ws.append(cells)
    ws.auto_filter.ref = f"A1:{get_column_letter(len(headers))}{row_number}"
    return digest.hexdigest()


def make_summary(counts, resolver, input_path, input_hash, sheet):
    labels = {
        "records": "\u603b\u8bb0\u5f55\u6570", "columns": "\u4e3b\u8868\u5217\u6570", "unique_pmids": "\u4e0d\u540cPMID\u503c\u6570\uff08\u542b\u65e0\u6548/\u7a7a\u503c\uff09",
        "changed_rows": "\u6709\u5b57\u6bb5\u53d8\u66f4\u7684\u8bb0\u5f55\u6570", "pending_rows": "\u542b\u5f85\u6838\u4e8b\u9879\u7684\u8bb0\u5f55\u6570",
        "hgvs_unconverted_rows": "\u672a\u5f62\u6210\u72ec\u7acb\u86cb\u767dHGVS\u5f62\u5f0f\u7684\u8bb0\u5f55\u6570\uff08\u542b\u7f3a\u5931/DNA/\u6784\u5efa\u4f53\uff09",
        "cell_unmatched_rows": "\u7ec6\u80de\u540d\u79f0\u7f3a\u5931/\u672a\u5b8c\u5168\u5339\u914d/\u6b67\u4e49\u8bb0\u5f55\u6570",
        "mechanism_pending_rows": "\u5b9a\u4f4d\u673a\u5236\u5f85\u6838\u6216\u6807\u7b7e\u51b2\u7a81\u8bb0\u5f55\u6570",
        "title_unverified_rows": "\u6807\u9898\u672a\u6838\u5b9a\u8bb0\u5f55\u6570", "title_unverified_pmids": "\u6807\u9898\u672a\u6838\u5b9aPMID\u6570",
        "main_duplicate_keys": "\u4e3b\u8868\u91cd\u590d\u53cc\u952e\u7ec4\u6570", "main_duplicate_key_rows": "\u4e3b\u8868\u91cd\u590d\u53cc\u952e\u8bb0\u5f55\u6570",
        "aux_conflict_rows": "\u8f85\u52a9\u8bc1\u636e\u53cc\u952e\u51b2\u7a81\u5f71\u54cd\u8bb0\u5f55\u6570", "aux_loaded_keys": "\u672c\u6b21\u4fdd\u7559\u7684\u8f85\u52a9\u53cc\u952e\u6570",
    }
    return [["\u7a0b\u5e8f/\u89c4\u5219\u7248\u672c", VERSION], ["\u8f93\u5165\u6587\u4ef6", str(input_path)], ["\u8f93\u5165SHA256", input_hash],
            ["\u8f93\u5165\u5de5\u4f5c\u8868", sheet], ["\u5904\u7406UTC\u65f6\u95f4", datetime.now(timezone.utc).isoformat()],
            ["\u6807\u9898\u6a21\u5f0f", resolver.mode], ["\u7edf\u8ba1\u53e3\u5f84", "\u5404\u5f85\u6838\u7c7b\u522b\u53ef\u91cd\u53e0\uff0c\u4e0d\u53ef\u76f8\u52a0\uff1b\u683c\u5f0f\u8f6c\u6362\u4e0d\u7b49\u4e8e\u79d1\u5b66\u6838\u9a8c"],
            *[[label, counts[key]] for key, label in labels.items()]]


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_rules(counts, resolver, sheet, missing_columns):
    return [
        ["\u7a0b\u5e8f\u7248\u672c", VERSION, "\u5355\u72ec\u672c\u6587\u4ef6\uff1b\u4ec5openpyxl\u4f9d\u8d56"],
        ["\u8f93\u5165\u5de5\u4f5c\u8868", sheet, "\u53ea\u5904\u7406\u9009\u4e2d\u8868\uff1b\u539f\u59cb\u5de5\u4f5c\u7c3f\u4e0d\u4fee\u6539"],
        ["\u8bb0\u5f55\u6570", counts["records"], "\u65e0\u5220\u9664\u3001\u65e0\u81ea\u52a8\u53bb\u91cd\uff1b\u539fExcel\u884c\u548cid\u53ef\u8ffd\u6eaf"],
        ["\u7f3a\u5c11\u7684\u6807\u51c6\u5217", "; ".join(missing_columns), "\u7f3a\u5931\u5217\u6309\u7a7a\u503c\u8865\u9f50\uff0c\u4e0d\u7f16\u9020\u4fe1\u606f"],
        ["\u7c7b\u578b\u89c4\u8303", "15\u7c7b\u9879\u76ee\u5206\u7c7b\uff1b\u4ec5\u72ec\u7acb\u86cb\u767d\u8868\u8fbe\u5f0f/\u4e25\u683c\u7b49\u4ef7\u522b\u540d\u8986\u76d6\u539f\u7c7b\u578b\u3002\u4e0a\u4e0b\u6587\u3001\u5bf9\u7167\u3001\u5426\u5b9a\u4e0e\u4e0d\u786e\u5b9a\u53e5\u4fdd\u7559\u6765\u6e90\u7c7b\u522b\u5e76\u5f85\u6838\u3002\u4e0d\u540c\u5c42\u7ea7\u4e0d\u80fd\u6df7\u5408\u4f5c\u53d8\u5f02\u9891\u7387\u7edf\u8ba1\u3002", "https://hgvs-nomenclature.org/stable/recommendations/checklist/"],
        ["HGVS", "\u53ea\u89c4\u8303\u65e0\u6b67\u4e49\u5c55\u793a\u5f62\u5f0f\uff1b\u672a\u6838\u5bf9\u5e8f\u5217\u3001\u7f16\u53f7\u3001\u8f6c\u5f55\u672c\u6216\u5b9e\u6d4b/\u9884\u6d4b\u540e\u679c\u3002\u591a\u53d8\u5f02\u548c\u5fc5\u8981\u4e0a\u4e0b\u6587\u4fdd\u7559\u539f\u6587\u3002", "https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/"],
        ["\u540c\u533a\u5ba4\u5b9a\u4f4d", "\u4f18\u5148\u5df2\u6709\u7ec6\u8282\u3001\u8868\u5185\u7ed3\u6784\u5316\u4e8b\u4ef6\uff1b\u539f\u59cb\u5b9a\u4f4d\u63cf\u8ff0\u4ec5\u4f5c\u5f85\u6838\u63d0\u793a\u3002\u6ca1\u6709\u4fe1\u606f\u5219\u6807\u5f85\u660e\u786e\uff0c\u4e0d\u5047\u8bbe\u6ede\u7559/\u805a\u96c6/\u6bd4\u4f8b\u53d8\u5316\u3002", "\u4e0d\u4f9d\u8d56\u5916\u90e8\u6e90JSON\uff1b\u65e0\u53d8\u5316\u4e0d\u7b49\u4e8e\u672a\u62a5\u544a"],
        ["\u7ec6\u80de\u540d\u79f0", "\u5185\u7f6e\u5e38\u89c1\u540d\u79f0\u7684\u4fdd\u5b88\u767d\u540d\u5355\uff0c\u5feb\u71672026-11-03\uff1b\u7cbe\u786e\u5b8c\u6574\u9879\u5339\u914d\uff0c\u4e0d\u5408\u5e76\u4e9a\u7cfb\u3002PC12/PC-12\u7b49\u6b67\u4e49\u4e0d\u81ea\u52a8\u6620\u5c04\u3002", "https://api.cellosaurus.org/api-fields"],
        ["\u6807\u9898", "\u6a21\u5f0f=" + resolver.mode + "\uff1b\u6309PMID\u5916\u90e8\u6838\u5b9a\u3002\u7f51\u7edc\u5931\u8d25/\u79bb\u7ebf\u672a\u547d\u4e2d\u7f13\u5b58\u5219\u539f\u503c\u4fdd\u7559\uff0c\u660e\u786e\u6807\u5f85\u6838\u3002", "https://europepmc.org/RestfulWebService"],
        ["\u9884\u5370\u672c\u6807\u9898", "\u901a\u8fc7PMC\u5b98\u65b9PMID\u2192DOI\u6620\u5c04\uff0c\u518d\u7cbe\u786e\u6838\u5bf9DOI\u539f\u7248\u672c\uff0c\u4e0d\u6539\u4e3a\u6b63\u5f0f\u51fa\u7248\u7248PMID\u3002", "https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api/"],
        ["\u7559\u7a7a\u5217", "; ".join(sorted(BLANK_COLUMNS)), "\u539f\u6709\u503c\u5728\u539f\u59cb\u6570\u636e\u53ca\u9010\u5b57\u6bb5\u5ba1\u8ba1\u4fdd\u7559"],
        ["\u75be\u75c5/\u57fa\u56e0/9\u5ba1\u9605\u5217", "\u4fdd\u7559\u8f93\u5165\u503c\uff0c\u53ea\u6807\u4e00\u81f4\u6027\u95ee\u9898\u3002\u6807\u7b7eNONE\u4e0e\u672a\u77e5\u4e0d\u540c\uff0c\u8f6c\u67d3\u4e0d\u81ea\u52a8\u7b49\u4e8e\u8fc7\u8868\u8fbe\uff1bNO/NOT_REPORTED/NOT_APPLICABLE\u4e0d\u6df7\u7528\u3002", "\u9879\u76ee\u5ba1\u9605\u5b57\u6bb5\u89c4\u5219"],
        ["\u8bc1\u636e\u5206\u7ea7", "\u9879\u76eeEvidence_Type\u4e0d\u7b49\u540c\u4e8eGO\u8bc1\u636e\u4ee3\u7801\uff1b\u4e0d\u80fd\u4ec5\u51ed\u65b9\u6cd5\u540d\u79f0\u5347\u7ea7\u8bc1\u636e\u6216\u8bad\u7ec3\u6807\u7b7e\u3002", "https://www.geneontology.org/docs/guide-go-evidence-codes/"],
        ["\u63a7\u5236\u5b57\u7b26\u4e0e\u516c\u5f0f", "\u4e0d\u517c\u5bb9XML\u7684\u63a7\u5236\u5b57\u7b26\u8f6c\u4e3a\u53ef\u89c1\\uXXXX\uff1b\u516c\u5f0f\u4f5c\u4e3a\u539f\u6587\u5b57\u7b26\u4e32\u4fdd\u5b58\uff0c\u4e0d\u6c42\u503c\u3002\u4e0d\u622a\u65ad\u8d85\u957f\u5355\u5143\u683c\u3002", "\u683c\u5f0f\u5904\u7406\uff0c\u4e0d\u4fee\u8865\u79d1\u5b66\u5185\u5bb9"],
        ["\u7f51\u7edc\u9519\u8bef", "\n".join(resolver.errors) or "\u65e0", "\u5931\u8d25\u4e0d\u4f2a\u88c5\u6210\u5b8c\u6210\u6838\u5b9a"],
        ["\u8f85\u52a9\u8fde\u63a5\u4e0e\u590d\u8dd1", "\u4e3b\u8868\u53ca\u8f85\u52a9\u8868\u53cc\u952e\u5747\u552f\u4e00\u624d\u8fde\u63a5\uff1b\u5df2\u6709\u7ec6\u8282\u91cd\u65b0\u68c0\u67e5\uff0c\u5e76\u4fdd\u7559\u53ef\u9760\u8fde\u63a5\u7684\u5386\u53f2\u5f85\u6838\u4e8b\u9879\u3002\u6807\u9898\u6309\u672c\u6b21\u6838\u5b9a\u91cd\u65b0\u8ba1\u7b97\u3002", "\u4e0d\u80fd\u4ee5\u91cd\u590d\u8fd0\u884c\u66ff\u4ee3\u4eba\u5de5\u89e3\u9664\u5f85\u6838"],
        ["\u6807\u9898\u7f13\u5b58", f"\u9010\u6279\u539f\u5b50\u4fdd\u5b58\uff1brefresh={resolver.refresh}\uff1bmax_age_days={resolver.max_age_days}", "\u79bb\u7ebf\u59cb\u7ec8\u4e0d\u8bf7\u6c42\u7f51\u7edc\uff1b\u8fc7\u671f/\u5237\u65b0\u5931\u8d25\u4fdd\u7559\u539f\u6807\u9898\u5e76\u6807\u5f85\u6838"],
        *[[k, v, "\u5904\u7406\u8ba1\u6570"] for k, v in sorted(counts.items()) if k.startswith("changed_")],
    ]


def run(args, resolver=None):
    from openpyxl import Workbook, load_workbook
    input_path = Path(args.input).expanduser().resolve(strict=True)
    output_path = Path(args.output).expanduser().resolve() if args.output else input_path.with_name(input_path.stem + "_\u89c4\u8303\u5316.xlsx")
    if input_path == output_path or (output_path.exists() and os.path.samefile(input_path, output_path)):
        raise ValueError("\u8f93\u51fa\u8def\u5f84\u4e0d\u80fd\u4e0e\u8f93\u5165\u76f8\u540c\uff0c\u539f\u59cbXLSX\u5fc5\u987b\u4fdd\u7559")
    if output_path.exists() and not args.force:
        raise ValueError(f"\u8f93\u51fa\u5df2\u5b58\u5728\uff1a{output_path}\uff1b\u66f4\u6362\u8f93\u51fa\u540d\u6216\u663e\u5f0f\u4f7f\u7528 --force")
    if input_path.suffix.lower() != ".xlsx" or output_path.suffix.lower() != ".xlsx":
        raise ValueError("\u8f93\u5165\u548c\u8f93\u51fa\u5fc5\u987b\u4e3a.xlsx\uff08\u4e0d\u652f\u6301\u5b8f\u5de5\u4f5c\u7c3f\u6216xls\uff09")
    cache_path = Path(args.cache).expanduser().resolve() if args.cache else None
    if cache_path and cache_path in {input_path, output_path}:
        raise ValueError("--cache \u4e0d\u5f97\u6307\u5411\u8f93\u5165/\u8f93\u51fa\u5de5\u4f5c\u7c3f")
    before_hash = file_sha256(input_path)
    wb = load_workbook(input_path, read_only=True, data_only=False)
    try:
        ws = select_sheet(wb, args.sheet)
        sheet_name = ws.title
        headers, source_rows = read_table(ws)
        wanted_keys = {(text(r.get("id")).strip(), text(r.get("PMID")).strip()) for _, r in source_rows if known(r.get("id"))}
        aux, conflicts = auxiliary_index(wb, sheet_name, wanted_keys)
    finally:
        wb.close()
    print(f"\u8bfb\u53d6 {sheet_name}: {len(source_rows)} \u6761\uff1b\u8f85\u52a9\u8bb0\u5f55 {len(aux)} \u6761", flush=True)
    mode = "cache" if args.offline else args.titles
    resolver = resolver or TitleResolver(mode, cache_path, args.timeout, args.retries, args.refresh_titles, args.cache_max_age_days)
    pmids = {text(r.get("PMID")).strip() for _, r in source_rows}
    resolver.resolve(pmids)
    unresolved = [p for p in pmids if resolver.mode == "keep" or not resolver.usable(p)]
    if args.strict_titles and unresolved:
        raise ValueError(f"\u8fd8\u6709{len(unresolved)}\u4e2aPMID\u672a\u6838\u5b9a\uff0c--strict-titles\u6a21\u5f0f\u4e0d\u5bfc\u51fa\uff1b\u793a\u4f8b\uff1a{unresolved[:10]}")
    out_headers, rows, changes, statuses, mappings, counts = normalize_rows(headers, source_rows, aux, conflicts, resolver)
    counts["title_unverified_pmids"] = len(unresolved)
    counts["aux_loaded_keys"] = len(aux)
    title_rows = []
    for p in sorted(pmids):
        e = resolver.entries.get(p, {}) if resolver.mode != "keep" and resolver.usable(p) else {}
        title_rows.append([p, e.get("title"), e.get("method", "\u672a\u5916\u90e8\u6838\u5b9a\uff0c\u4fdd\u7559\u9010\u6761\u539f\u6807\u9898"), e.get("url"), e.get("binding_url"), e.get("retrieved_at")])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=output_path.stem + ".", suffix=".xlsx", dir=output_path.parent)
    os.close(handle)
    expected = {}
    output_wb = Workbook(write_only=True)
    try:
        expected["\u5904\u7406\u6458\u8981"] = write_sheet(output_wb, "\u5904\u7406\u6458\u8981", ["\u6307\u6807", "\u503c"], make_summary(counts, resolver, input_path, before_hash, sheet_name), [75, 115])
        expected["\u89c4\u8303\u5316\u6570\u636e"] = write_sheet(output_wb, "\u89c4\u8303\u5316\u6570\u636e", out_headers, ([r[h] for h in out_headers] for r in rows))
        expected["\u539f\u59cb\u6570\u636e"] = write_sheet(output_wb, "\u539f\u59cb\u6570\u636e", ["\u6765\u6e90Excel\u884c", *headers], ([i, *[r[h] for h in headers]] for i, r in source_rows))
        expected["\u9010\u5b57\u6bb5\u5ba1\u8ba1"] = write_sheet(output_wb, "\u9010\u5b57\u6bb5\u5ba1\u8ba1", ["\u6765\u6e90Excel\u884c", "id", "PMID", "\u5b57\u6bb5", "\u539f\u503c", "\u89c4\u8303\u503c", "\u4f9d\u636e"], changes, [14, 14, 14, 40, 80, 80, 80])
        expected["\u5904\u7406\u72b6\u6001"] = write_sheet(output_wb, "\u5904\u7406\u72b6\u6001", ["\u6765\u6e90Excel\u884c", "id", "PMID", "HGVS\u5904\u7406\u72b6\u6001", "\u5b9a\u4f4d\u7ec6\u8282\u4f9d\u636e", "Cellosaurus\u7f16\u53f7", "\u6807\u9898\u72b6\u6001", "\u5f85\u6838\u4e8b\u9879"], statuses, [14, 14, 14, 70, 80, 30, 60, 100])
        expected["\u7c7b\u578b\u6620\u5c04"] = write_sheet(output_wb, "\u7c7b\u578b\u6620\u5c04", ["\u539f\u7c7b\u578b", "\u89c4\u8303\u7c7b\u578b", "\u8bb0\u5f55\u6570", "\u4f9d\u636e"], [[a, b, n, why] for (a, b, why), n in sorted(mappings.items())], [50, 30, 15, 80])
        expected["\u6807\u9898\u6838\u5b9a"] = write_sheet(output_wb, "\u6807\u9898\u6838\u5b9a", ["PMID", "\u6838\u5b9a\u6807\u9898", "\u65b9\u5f0f", "\u6765\u6e90", "PMID-DOI\u6838\u5b9a\u6765\u6e90", "\u67e5\u8be2UTC\u65f6\u95f4"], title_rows, [15, 100, 65, 75, 75, 30])
        expected["\u89c4\u5219\u8bf4\u660e"] = write_sheet(output_wb, "\u89c4\u5219\u8bf4\u660e", ["\u9879\u76ee", "\u89c4\u5219/\u503c", "\u6765\u6e90/\u5907\u6ce8"], make_rules(counts, resolver, sheet_name, [h for h in BASE_COLUMNS if h not in headers]), [35, 115, 85])
        print("\u5199\u5165\u5e76\u9010\u5355\u5143\u683c\u56de\u8bfb\u6821\u9a8c\u2026\u2026", flush=True)
        output_wb.save(temporary)
        check = load_workbook(temporary, read_only=True, data_only=False)
        try:
            for name, digest in expected.items():
                if data_digest(check[name].values) != digest:
                    raise ValueError(f"\u4fdd\u5b58\u540e\u5185\u5bb9\u6821\u9a8c\u5931\u8d25\uff1a{name}")
        finally:
            check.close()
        if file_sha256(input_path) != before_hash:
            raise ValueError("\u5904\u7406\u671f\u95f4\u8f93\u5165\u6587\u4ef6\u53d1\u751f\u53d8\u5316\uff0c\u505c\u6b62\u4ea4\u4ed8\uff0c\u8bf7\u91cd\u65b0\u8fd0\u884c")
        if output_path.exists() and not args.force:
            raise ValueError("\u5904\u7406\u671f\u95f4\u8f93\u51fa\u8def\u5f84\u88ab\u5360\u7528\uff0c\u8bf7\u6362\u8f93\u51fa\u540d")
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    print(f"\u5b8c\u6210\uff1a{output_path}\n\u8bb0\u5f55\u6570\uff1a{len(rows)}\uff1b\u5217\u6570\uff1a{len(out_headers)}\uff1b\u672a\u6838\u5b9a\u6807\u9898PMID\uff1a{len(unresolved)}", flush=True)
    return {"output": str(output_path), "rows": len(rows), "columns": len(out_headers), "unresolved_titles": len(unresolved), "counts": dict(counts)}


def self_test():
    """\u72ec\u7acb\u81ea\u68c0\uff1a\u7eaf\u51fd\u6570\u3001\u7f51\u7edc\u54cd\u5e94\u8fb9\u754c\u548c\u771f\u6b63\u7684XLSX\u8f93\u5165/\u8f93\u51fa\u5747\u6d4b\u8bd5\uff0c\u4e0d\u8054\u7f51\u3002"""
    import unittest
    from openpyxl import Workbook, load_workbook

    class Tests(unittest.TestCase):
        def test_protein_cases(self):
            self.assertEqual(len(AA), 20)
            self.assertEqual(AA['P'], 'Pro')
            self.assertEqual(AA['V'], 'Val')
            cases = [("V107R", "Missense", "p.Val107Arg", "Missense"),
                     ("F508del (p.Phe508del)", "Missense", "p.Phe508del", "Deletion"),
                     ("p.Phe508del (F508del)", "Missense", "p.Phe508del", "Deletion"),
                     ("p.Arg4218Ter", "Missense", "p.Arg4218Ter", "Nonsense"),
                     ("p.Arg97ProfsTer23", "Missense", "p.Arg97ProfsTer23", "Frameshift"),
                     ("p.(R123S)", "Missense", "p.(Arg123Ser)", "Missense"),
                     ("p.(Arg123Ser)", "Missense", "p.(Arg123Ser)", "Missense"),
                     ("p.Lys237_Val238del", "Missense", "p.Lys237_Val238del", "Deletion"),
                     ("p.Arg76_Cys77delinsSerTrp", "Missense", "p.Arg76_Cys77delinsSerTrp", "Delins"),
                     ("R123X", "Other", "R123X", "Other"),
                     ("R123X", "Nonsense", "p.Arg123Ter", "Nonsense"),
                     ("c.123+1G>A", "Splice-site", "c.123+1G>A", "Splice-site")]
            for raw, old, value, typ in cases:
                with self.subTest(raw=raw):
                    self.assertEqual(normalize_mutation(raw, old)[0], value)
                    self.assertEqual(normalize_type(raw, old)[0], typ)

        def test_ambiguity_and_negation(self):
            self.assertEqual(normalize_type("p.=", "Missense")[0], "Unresolved")
            self.assertEqual(normalize_type("p.Arg123Arg", "Missense")[0], "Unresolved")
            self.assertEqual(normalize_type("p.Arg123=", "Other")[0], "Unresolved")
            self.assertEqual(normalize_type("p.Met1Val", "Missense")[0], "Unresolved")
            self.assertEqual(normalize_type("in-frame deletion; no frameshift", "Deletion")[0], "Deletion")
            self.assertEqual(normalize_type("Wild-type Gish-GFP; no mutation", "None")[0], "Wild-type")
            self.assertEqual(normalize_type("siRNA knockdown", "Knockdown")[0], "Non-sequence perturbation")
            value = "p.Arg12Gly; p.Leu20del (different constructs)"
            self.assertEqual(normalize_mutation(value, "Missense")[0], value)
            self.assertEqual(normalize_type(value, "Missense")[0], "Missense")
            self.assertIn("\u5f85\u6838", normalize_type(value, "Missense")[1])
            self.assertIsNone(short_protein("R0W"))

        def test_context_does_not_reclassify(self):
            for raw in ["Deletion of exon 4; compared with p.Arg123Trp control",
                        "not a frameshift; p.Arg123Trp comparator", "no evidence of frameshift",
                        "frameshift cannot be excluded", "possible frameshift", "without any frameshift",
                        "non-frameshift deletion", "p.Arg123Trp versus p.Leu20del"]:
                with self.subTest(raw=raw):
                    self.assertEqual(normalize_type(raw, "Deletion")[0], "Deletion")
                    self.assertIn("\u5f85\u6838", normalize_type(raw, "Deletion")[1])
            self.assertEqual(normalize_type("frameshift mutation", "Deletion")[0], "Frameshift")

        def test_localization_conflicts_and_second_pass(self):
            base = {"Mislocalization Phenotype\uff08before\uff09": "Nucleus", "Mislocalization Phenotype\uff08after\uff09": "Nucleus"}
            for outcome, event, extent in [("EXPERIMENTAL_NO_CHANGE", "aggregation", ""),
                                           ("EXPERIMENTAL_MISLOCALIZATION", "no_change", ""),
                                           ("", "retention", "NO_CHANGE"),
                                           ("", "translocation", ""),
                                           ("EXPERIMENTAL_NO_CHANGE", "", "QUANTITATIVE")]:
                row = dict(base, Outcome_Provenance=outcome, Redistribution_Extent=extent)
                detail, _, flags = localization_detail(row, {"localization_event_type": event})
                self.assertIn("\u51b2\u7a81", detail)
                self.assertTrue(flags)
                rerun = localization_detail(dict(row, **{DETAIL: detail}), {})
                self.assertEqual(rerun[0], detail)
                self.assertTrue(rerun[2])
            changed = dict(base, **{"Mislocalization Phenotype\uff08after\uff09": "Cytoplasm", "Outcome_Provenance": "EXPERIMENTAL_NO_CHANGE"})
            self.assertIn("\u51b2\u7a81", localization_detail(changed, {})[0])
            for detail in ["\u673a\u5236\u5f85\u660e\u786e", "\u540c\u533a\u5ba4\u53d8\u5316\uff0c\u673a\u5236\u5f85\u660e\u786e", "\u518d\u5206\u5e03\uff08\u5177\u4f53\u5f62\u5f0f\u5f85\u6838\uff09", "\u539f\u59cb\u5b9a\u4f4d\u63cf\u8ff0\u63d0\u793a\u805a\u96c6\uff08\u5f85\u590d\u6838\uff09"]:
                self.assertTrue(localization_detail(dict(base, **{DETAIL: detail}), {})[2])
            self.assertEqual(localization_detail(base, {"localization_event_type": "retention"})[0], "\u6ede\u7559")
            conflict = dict(base, **{DETAIL: "\u805a\u96c6", "Outcome_Provenance": "EXPERIMENTAL_NO_CHANGE"})
            self.assertTrue(localization_detail(conflict, {})[2])

        def test_structured_tags(self):
            for value in ["kind=FLUORESCENT_PROTEIN; position = NONE; role=DETECTION_REPORTER",
                          " KIND = epitope_affinity \uff1b POSITION = none "]:
                row = {"Tag_Strategy": value}
                self.assertIn("\u6709\u6807\u7b7e\u4f46\u4f4d\u7f6eNONE\uff0c\u9700\u67e5\u539f\u6587\uff1b\u672a\u81ea\u52a8\u6539\u5199", review_flags(row))
                self.assertEqual(row["Tag_Strategy"], value)
            self.assertTrue(review_flags({"Tag_Strategy": "position=N; position=C"}))

        def test_main_duplicate_keys_disable_auxiliary(self):
            base = {"id": 1, "PMID": "123", "Mutation Type": "Missense", "Mislocalization Phenotype\uff08before\uff09": "Nucleus", "Mislocalization Phenotype\uff08after\uff09": "Nucleus"}
            source = [(2, dict(base, **{"Mutation Information": "R123W"})), (3, dict(base, **{"Mutation Information": "V25A"}))]
            result = normalize_rows(list(source[0][1]), source, {("1", "123"): {"localization_event_type": "retention"}}, set(), TitleResolver("cache"))
            self.assertTrue(all(r[DETAIL] == "\u540c\u533a\u5ba4\u53d8\u5316\uff0c\u673a\u5236\u5f85\u660e\u786e" for r in result[1]))
            self.assertTrue(all("\u4e3b\u8868\u53cc\u952e\u4e0d\u552f\u4e00" in r[-1] for r in result[3]))
            self.assertEqual(result[5]["main_duplicate_keys"], 1)
            self.assertEqual(result[5]["main_duplicate_key_rows"], 2)

        def test_cache_checkpoint_resume_expiry_refresh(self):
            from unittest.mock import patch
            class Interrupt(TitleResolver):
                def get_json(self, base, params):
                    ids = re.findall(r"EXT_ID:(\d+)", params.get("query", ""))
                    if "81" in ids:
                        raise KeyboardInterrupt()
                    return {"resultList": {"result": [{"source": "MED", "id": p, "title": "Title " + p} for p in ids]}}
            class NoNetwork(TitleResolver):
                def get_json(self, base, params):
                    raise AssertionError("\u4e0d\u5e94\u8bf7\u6c42\u7f51\u7edc")
            with tempfile.TemporaryDirectory() as td, patch("time.sleep"):
                cache = Path(td) / "titles.json"
                resolver = Interrupt(cache_path=cache)
                with self.assertRaises(KeyboardInterrupt):
                    resolver.resolve({str(i) for i in range(1, 82)})
                self.assertEqual(len(json.loads(cache.read_text(encoding="utf-8"))["titles"]), 80)
                offline = NoNetwork("cache", cache)
                offline.resolve({"1", "81"})
                self.assertTrue(offline.usable("1"))
                self.assertFalse(offline.usable("81"))
                resume = NoNetwork("online", cache)
                resume.resolve({"1", "80"})  # \u5df2\u5b8c\u6210\u7684\u6279\u6b21\u4e0d\u518d\u8054\u7f51\u3002
                expired = NoNetwork("cache", cache, max_age_days=1)
                expired.entries["1"]["retrieved_at"] = "2000-01-01T00:00:00+00:00"
                self.assertFalse(expired.usable("1"))
                refresh = TitleResolver("online", cache, refresh=True)
                self.assertFalse(refresh.usable("1"))
                refresh.remember("1", {"title": "New title"}, "mock", "https://example.org/1")
                self.assertTrue(refresh.usable("1"))
                self.assertTrue(NoNetwork("cache", cache, refresh=True).usable("1"))

        def test_auxiliary_stream_filter_and_prior_flags(self):
            wb = Workbook()
            ws = wb.create_sheet("\u590d\u6838\u7d22\u5f15")
            ws.append(["id", "PMID", "localization_event_type"])
            ws.append([1, 123, "retention"])
            ws.append([2, 123, "aggregation"])
            ws.append([2, 123, "retention"])
            prior = wb.create_sheet("\u5904\u7406\u72b6\u6001")
            prior.append(["id", "PMID", "\u5f85\u6838\u4e8b\u9879", "\u5b9a\u4f4d\u7ec6\u8282\u4f9d\u636e"])
            prior.append([1, 123, "\u5386\u53f2\u8bc1\u636e\u5f85\u6838; \u6807\u9898\u5f85\u5916\u90e8\u6838\u5b9a", "\u6e90\u7ed3\u6784\u5316\u4e8b\u4ef6"])
            index, conflicts = auxiliary_index(wb, "Sheet", {("1", "123")})
            self.assertEqual(set(index), {("1", "123")})
            self.assertFalse(conflicts)
            row = {"id": 1, "PMID": "123", DETAIL: "\u6ede\u7559"}
            resolver = TitleResolver("cache")
            resolver.remember("123", {"title": "Verified"}, "mock", "https://example.org/123")
            result = normalize_rows(list(row), [(2, row)], index, conflicts, resolver)
            self.assertIn("\u5386\u53f2\u8bc1\u636e\u5f85\u6838", result[3][0][-1])
            self.assertNotIn("\u6807\u9898\u5f85\u5916\u90e8\u6838\u5b9a", result[3][0][-1])
            self.assertEqual(result[3][0][4], "\u6e90\u7ed3\u6784\u5316\u4e8b\u4ef6")

        def test_cells(self):
            self.assertEqual(normalize_cell("HEK 293; HEK-293T")[0], "HEK293; HEK293T")
            self.assertEqual(normalize_cell("HEK;293")[0], "HEK293")
            self.assertEqual(normalize_cell("HeLa (Kyoto)")[0], "HeLa Kyoto")
            self.assertEqual(normalize_cell("HEK293T/17")[0], "HEK293T/17")
            self.assertEqual(normalize_cell("PC12")[1], [])
            self.assertEqual(normalize_cell("HEK293 expressing mutant ABC")[0], "HEK293 expressing mutant ABC")

        def test_no_invented_mechanism(self):
            row = {"Mislocalization Phenotype\uff08before\uff09": "Nucleus", "Mislocalization Phenotype\uff08after\uff09": "Nucleus", "Outcome_Provenance": "EXPERIMENTAL_MISLOCALIZATION"}
            self.assertEqual(localization_detail(row, {})[0], "\u540c\u533a\u5ba4\u53d8\u5316\uff0c\u673a\u5236\u5f85\u660e\u786e")
            self.assertEqual(localization_detail(row, {"localization_event_type": "retention"})[0], "\u6ede\u7559")
            self.assertIn("\u51b2\u7a81", localization_detail(row, {"localization_event_type": "no_change"})[0])
            self.assertNotIn("\u63d0\u793a", localization_detail(row, {"mutated_localization_raw": "no aggregation"})[0])

        def test_missing_id_and_auxiliary_collision(self):
            row = {"PMID": "123", "Mutation Information": "R123W", "Mutation Type": "Missense", "Mislocalization Phenotype": "Nucleus>Nucleus", "Mislocalization Phenotype\uff08before\uff09": "Nucleus", "Mislocalization Phenotype\uff08after\uff09": "Nucleus"}
            result = normalize_rows(list(row), [(2, row)], {("xlsx-row-2", "123"): {"localization_event_type": "retention"}}, set(), TitleResolver("cache"))
            self.assertEqual(result[1][0][DETAIL], "\u540c\u533a\u5ba4\u53d8\u5316\uff0c\u673a\u5236\u5f85\u660e\u786e")
            self.assertEqual(result[1][0]["id"], "xlsx-row-2")
            wb = Workbook()
            ws = wb.create_sheet("\u590d\u6838\u7d22\u5f15")
            ws.append(["id", "PMID", "localization_event_type"])
            ws.append([1, 123, "retention"])
            ws.append([1, 123, "aggregation"])
            index, conflicts = auxiliary_index(wb, "Sheet")
            self.assertIn(("1", "123"), conflicts)
            self.assertNotIn(("1", "123"), index)

        def test_title_wrong_id_and_network_failure(self):
            class Mock(TitleResolver):
                def get_json(self, base, params):
                    if "idconv" in base:
                        return {"records": []}
                    return {"resultList": {"result": [{"source": "MED", "id": "123", "pmid": "123", "title": "Title"}, {"source": "MED", "id": "999", "pmid": "456", "title": "Wrong binding"}]}}
            r = Mock()
            r.resolve({"123", "456"})
            self.assertTrue(r.usable("123"))
            self.assertFalse(r.usable("456"))
            class Down(TitleResolver):
                def get_json(self, base, params):
                    raise RuntimeError("offline")
            r = Down()
            r.resolve({"123"})
            self.assertFalse(r.usable("123"))
            self.assertTrue(r.errors)

        def test_workbook_end_to_end(self):
            with tempfile.TemporaryDirectory(prefix="xlsx_norm_test_") as td:
                src, dst = Path(td) / "input.xlsx", Path(td) / "out.xlsx"
                wb = Workbook()
                sh = wb.active
                sh.title = "Data"
                sh.append(BASE_COLUMNS)
                row = {"id": 10, "PMID": 123, "Mutation Type": "Missense", "Mutation Information": "F508del", "Gene Name": "CFTR", "Mislocalization Phenotype": "ER>ER", "Mislocalization Phenotype\uff08before\uff09": "ER", "Mislocalization Phenotype\uff08after\uff09": "ER", "Cell Line": "HEK-293", "Disease Information": "keep me", "Title": "original", "Accession_new": "old accession", "Outcome_Provenance": "EXPERIMENTAL_MISLOCALIZATION"}
                sh.append([row.get(h) for h in BASE_COLUMNS])
                pending = dict(row, id=20, **{"Mutation Type": "Deletion", "Mutation Information": "no evidence of frameshift", "Tag_Strategy": "kind=FLUORESCENT_PROTEIN; position = NONE"})
                sh.append([pending.get(h) for h in BASE_COLUMNS])
                aux = wb.create_sheet("\u590d\u6838\u7d22\u5f15")
                aux.append(["id", "PMID", "localization_event_type"])
                aux.append([10, 123, "retention"])
                wb.save(src)
                before = src.read_bytes()
                args = parser().parse_args([str(src), "-o", str(dst), "--offline"])
                run(args)
                check = load_workbook(dst, data_only=False)
                values = list(check["\u89c4\u8303\u5316\u6570\u636e"].values)
                result = dict(zip(values[0], values[1]))
                self.assertEqual(result[DETAIL], "\u6ede\u7559")
                self.assertEqual(result["Mutation Type"], "Deletion")
                self.assertEqual(result["Mutation Information"], "p.Phe508del")
                self.assertEqual(result["Cell Line"], "HEK293")
                self.assertEqual(result["Disease Information"], "keep me")
                self.assertTrue(all(result[h] is None for h in BLANK_COLUMNS))
                self.assertEqual(src.read_bytes(), before)
                first_status = list(check["\u5904\u7406\u72b6\u6001"].values)
                summary = dict(list(check["\u5904\u7406\u6458\u8981"].values)[1:])
                self.assertEqual(summary["\u603b\u8bb0\u5f55\u6570"], 2)
                self.assertEqual(summary["\u8f93\u5165SHA256"], file_sha256(src))
                check.close()
                rerun_path = Path(td) / "rerun.xlsx"
                run(parser().parse_args([str(dst), "-o", str(rerun_path), "--offline"]))
                rerun = load_workbook(rerun_path, read_only=True)
                self.assertEqual(list(rerun["\u89c4\u8303\u5316\u6570\u636e"].values), values)
                for first, second in zip(first_status[1:], list(rerun["\u5904\u7406\u72b6\u6001"].values)[1:]):
                    self.assertLessEqual(set(first[-1].split("; ")), set(second[-1].split("; ")))
                rerun.close()
                with self.assertRaises(ValueError):
                    run(args)  # \u62d2\u7edd\u9759\u9ed8\u8986\u76d6\u3002
                args.output = str(src)
                with self.assertRaises(ValueError):
                    run(args)

        def test_literal_formula_and_control(self):
            self.assertEqual(excel_value("a\x01b"), "a\\u0001b")
            with self.assertRaises(ValueError):
                excel_value("x" * 32768)
            with tempfile.TemporaryDirectory() as td:
                path = Path(td) / "literal.xlsx"
                wb = Workbook(write_only=True)
                expected = write_sheet(wb, "Data", ["value"], [["=1+1"], ["a\x01b"]])
                wb.save(path)
                check = load_workbook(path, read_only=True, data_only=False)
                self.assertEqual(data_digest(check["Data"].values), expected)
                self.assertEqual(list(check["Data"].rows)[1][0].data_type, "s")
                check.close()

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    return unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="?", help="\u8f93\u5165XLSX\u5b8c\u6574\u8def\u5f84")
    p.add_argument("-o", "--output", help="\u8f93\u51faXLSX\uff1b\u9ed8\u8ba4\u8f93\u5165\u540d\u52a0_\u89c4\u8303\u5316")
    p.add_argument("--sheet", help="\u6570\u636e\u5de5\u4f5c\u8868\u540d\u79f0\uff1b\u7b2c\u4e00\u884c\u4e3a\u8868\u5934")
    p.add_argument("--titles", choices=["online", "cache", "keep"], default="online", help="\u6807\u9898\u6a21\u5f0f\uff0c\u9ed8\u8ba4online\uff1bcache\u53ea\u8bfb\u7f13\u5b58\uff0ckeep\u4e0d\u6539\u6807\u9898")
    p.add_argument("--offline", action="store_true", help="\u4e0d\u8054\u7f51\uff0c\u7b49\u6548--titles cache\uff1b\u65e0\u7f13\u5b58\u4fdd\u7559\u539f\u6807\u9898")
    p.add_argument("--cache", help="\u53ef\u9009\u6807\u9898\u7f13\u5b58JSON\uff1b\u4e0d\u6307\u5b9a\u5219\u53ea\u4f7f\u7528\u5185\u5b58\uff0c\u65e0\u5916\u90e8\u6587\u4ef6\u4f9d\u8d56")
    p.add_argument("--refresh-titles", action="store_true", help="\u5728\u7ebf\u65f6\u91cd\u65b0\u6838\u5b9a\u5df2\u7f13\u5b58\u6807\u9898\uff1b\u79bb\u7ebf\u65f6\u65e0\u6548")
    p.add_argument("--cache-max-age-days", type=float, help="\u7f13\u5b58\u6709\u6548\u671f\u5929\u6570\uff1b\u9ed8\u8ba4\u4e0d\u9650\uff0c\u8fc7\u671f\u503c\u4e0d\u89c6\u4e3a\u5df2\u6838\u5b9a\uff1b\u79bb\u7ebf\u4e5f\u6267\u884c\u671f\u9650\u68c0\u67e5")
    p.add_argument("--strict-titles", action="store_true", help="\u4efb\u4e00PMID\u6807\u9898\u672a\u5916\u90e8\u6838\u5b9a\u5219\u4e0d\u5bfc\u51fa")
    p.add_argument("--timeout", type=float, default=20, help="\u5355\u6b21\u7f51\u7edc\u8bf7\u6c42\u8d85\u65f6\u79d2\u6570\uff0c\u9ed8\u8ba420")
    p.add_argument("--retries", type=int, default=1, choices=range(0, 4), help="\u5355\u8bf7\u6c42\u91cd\u8bd5\u6b21\u6570\uff0c\u9ed8\u8ba41\uff0c\u6700\u59273")
    p.add_argument("--force", action="store_true", help="\u5141\u8bb8\u66ff\u6362\u5df2\u6709\u8f93\u51fa\u6587\u4ef6\uff0c\u59cb\u7ec8\u7981\u6b62\u8986\u76d6\u8f93\u5165")
    p.add_argument("--self-test", action="store_true", help="\u8fd0\u884c\u5185\u7f6e\u6d4b\u8bd5\uff0c\u4e0d\u8054\u7f51\u4e14\u4e0d\u9700\u8f93\u5165\u6587\u4ef6")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def main():
    args = parser().parse_args()
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        print("\u7f3a\u5c11\u4f9d\u8d56\uff0c\u8bf7\u8fd0\u884c\uff1apython -m pip install openpyxl", file=sys.stderr)
        return 2
    if args.self_test:
        return 0 if self_test() else 2
    if not args.input:
        parser().error("\u8bf7\u6307\u5b9a\u8f93\u5165XLSX\uff0c\u6216\u4f7f\u7528 --self-test")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser().error("--timeout \u5fc5\u987b\u4e3a\u6b63\u6570")
    if args.cache_max_age_days is not None and (not math.isfinite(args.cache_max_age_days) or args.cache_max_age_days <= 0):
        parser().error("--cache-max-age-days \u5fc5\u987b\u4e3a\u6709\u9650\u6b63\u6570")
    try:
        run(args)
        return 0
    except KeyboardInterrupt:
        print("\u5df2\u4e2d\u65ad\uff1b\u8f93\u5165\u6587\u4ef6\u672a\u4fee\u6539\u3002", file=sys.stderr)
        return 130
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"\u5931\u8d25\uff1a{exc}", file=sys.stderr)
        return 2



CORE_SOURCE_SHA256 = '57a8397ad737c429bde8d003d6ca2fe6d4dbe419f09893b2773bcca2411deba4'

# Step 05 JSON adapter. Core normalization rules above are frozen at 1.1.0.
ROOT = Path(__file__).resolve().parents[2]

def is_human_target(record):
    """Use the target protein species only, never the expression host."""
    raw_value = re.sub(r'\s+', ' ', text(record.get('target_species')).strip())
    value = raw_value.casefold()
    if value in {'human', 'homo sapiens', 'homo sapiens (human)', 'human (homo sapiens)', '9606'}:
        return True
    # Accept an explicit taxonomy label with a benign qualifier, not any phrase
    # containing "human" (which would include human viruses and host-cell names).
    match = re.fullmatch(r'(?:homo sapiens|human)\s*\(([^()]+)\)', value)
    if not match:
        return False
    qualifier = match[1].strip()
    if re.search(r'virus|viral|hiv|hsv|hcmv|papilloma|mouse|murine|rat\b|bovine|yeast|'
                 r'hamster|monkey|zebrafish|canine|feline|dog\b|cat\b|non[- ]?human|chimer|mixed|unknown|inferred|likely|presumed|host|cell|\bor\b|\band\b|[;/]', qualifier):
        return False
    if re.fullmatch(r'(?:homo sapiens|human|(?:tax(?:on|onomy)?(?:\s*id)?\s*[:=]?\s*)?9606)', qualifier):
        return True
    raw_qualifier = raw_value[raw_value.rfind('(')+1:-1].strip()
    return bool(re.fullmatch(r'[A-Z][A-Z0-9-]{1,24}',raw_qualifier))


def explicit_conflict_reasons(record, normalization_flags):
    """Exclude positive conflict flags, not unresolved identity or lookup misses."""
    reasons = []
    for flag in text(normalization_flags).split('; '):
        if '\u51b2\u7a81' in flag:
            reasons.append('normalization_conflict: ' + flag)
    if text(record.get('relationship_verification')).upper() == 'CONFLICT':
        reasons.append('protein_relationship_conflict')
    if 'evidence_method_conflict' in text(record.get('error_type')).split(';'):
        reasons.append('evidence_method_conflict')
    return list(dict.fromkeys(reasons))


def uniprot_query_name(value):
    name = text(value).strip()
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9 .-]{0,60}', name):
        return None
    if name.casefold() in UNKNOWN or re.search(
            r'\b(?:mutant|mutated|wild|fusion|tagged|variant|construct|deletion|truncated)\b|'
            r'(?:^|[- ])(?:GFP|EGFP|YFP|ECFP|FLAG|HA|MYC)(?:[- ]|$)|'
            r'(?<![A-Za-z])[A-Z][0-9]{2,}[A-Z](?![A-Za-z])', name, re.I):
        return None
    return name


def strict_quality_decision(record, source_text, base_mapping='SKIPPED'):
    """Route unsupported claims to review, never describe them as false data."""
    attribution = text(record.get('source_attribution')).casefold()
    quote = text(record.get('evidence_quote'))
    secondary = (record.get('Outcome_Provenance') == 'SECONDARY_ONLY' or
                 record.get('source_attribution_class') == 'SECONDARY_ONLY' or
                 bool(re.search(r'not this study.{0,25}(?:experiment|primary)|citing |attributed to reference|'
                                r'background only|secondary evidence', attribution)))
    if secondary:
        return 'OUT_OF_SCOPE', ['secondary_source_not_primary_experiment']
    if base_mapping == 'DISTINCT_BASE_PROTEINS' or record.get('variant_observed_relationship') == 'DOWNSTREAM_EFFECT':
        return 'OUT_OF_SCOPE', ['observed_protein_is_not_mutated_protein']
    reasons = []
    if record.get('Evidence_Type') not in {'SPATIAL_DIRECT', 'BIOCHEMICAL_LOCALIZATION'}:
        reasons.append('indirect_or_unconfirmed_localization_evidence')
    if record.get('outcome_class') not in {'mislocalization', 'normal_experimentally_confirmed'}:
        reasons.append('outcome_not_confirmed')
    if record.get('variant_observed_relationship') != 'SAME_PROTEIN':
        reasons.append('mutant_observed_identity_unresolved')
    if record.get('relationship_verification') != 'SOURCE_IDENTITY_CONFIRMED':
        reasons.append('mutant_identity_not_source_confirmed')
    if not source_text:
        reasons.append('source_markdown_missing')
    else:
        norm = lambda value: re.sub(r'\s+', ' ', text(value)).strip()
        if not quote or norm(quote) not in norm(source_text):
            reasons.append('observation_quote_not_matched_to_source')
    if not re.search(r'\b(?:this|present|current) study\b|our (?:study|experiments|results)', attribution):
        reasons.append('primary_study_attribution_unconfirmed')
    assessment = record.get('evidence_assessment') or {}
    for key in ['gene_name', 'mutation_information', 'normal_localization',
                'mutated_localization', 'mislocalization_phenotype', 'outcome_class']:
        if assessment.get(key) != 'SUPPORTED':
            reasons.append('core_evidence_pending:' + key)
    return ('NEEDS_REVIEW', reasons) if reasons else ('HIGH_CONFIDENCE_CANDIDATE', [])


class HumanUniProtLookup:
    """Bounded exact-name lookup; API failure is distinct from no matching entry."""
    SCHEMA = 'human-uniprot-name-cache-v1'
    URL = 'https://rest.uniprot.org/uniprotkb/search'

    def __init__(self, cache_path, session=None, timeout=15):
        import requests
        self.session = session or requests.Session()
        self.cache_path, self.timeout = cache_path, timeout
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.entries, self.last_request, self.failure_streak = {}, 0.0, 0
        self.requests, self.cache_hits = 0, 0
        if cache_path.exists():
            content = json.loads(cache_path.read_text(encoding='utf-8'))
            if content.get('schema') != self.SCHEMA:
                raise ValueError('Unsupported UniProt cache schema.')
            self.entries = content['entries']

    @staticmethod
    def exact_names(entry):
        names = []
        for gene in entry.get('genes', []):
            names.append(gene.get('geneName', {}).get('value', ''))
            for field in ('synonyms', 'orderedLocusNames', 'orfNames'):
                names.extend(n.get('value', '') for n in gene.get(field, []))
        description = entry.get('proteinDescription', {})
        for section in [description.get('recommendedName', {}), *description.get('alternativeNames', [])]:
            names.append(section.get('fullName', {}).get('value', ''))
            names.extend(n.get('value', '') for n in section.get('shortNames', []))
        return names

    def lookup(self, value):
        name = uniprot_query_name(value)
        if not name:
            return {'status': 'SKIPPED_NAME', 'matches': []}
        key = name.casefold()
        cached = self.entries.get(key)
        # Cached misses are explicitly scoped to reviewed entries at lookup time.
        if cached and cached['status'] in {'UNIQUE', 'AMBIGUOUS', 'NOT_FOUND_REVIEWED', 'TRUNCATED'}:
            self.cache_hits += 1
            return cached
        if self.failure_streak >= 3:
            return {'status': 'DEFERRED_AFTER_REQUEST_FAILURES', 'matches': []}
        query = f'organism_id:9606 AND reviewed:true AND (gene_exact:"{name}" OR protein_name:"{name}")'
        result = {'status': 'REQUEST_FAILED', 'name': name, 'taxon_id': 9606,
                  'query': query, 'queried_at': datetime.now(timezone.utc).isoformat(), 'matches': []}
        try:
            time.sleep(max(0, 1.0 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            self.requests += 1
            response = self.session.get(self.URL, params={
                'query': query, 'format': 'json', 'size': 100,
                'fields': 'accession,gene_names,protein_name,organism_id,organism_name'},
                timeout=(5, self.timeout))
            response.raise_for_status()
            payload = response.json()
            entries = payload.get('results')
            if not isinstance(entries, list):
                raise ValueError('UniProt response has no results array.')
            result['response_url'] = response.url
            result['database_release'] = response.headers.get('X-UniProt-Release', '')
            for entry in entries:
                aliases = self.exact_names(entry)
                if entry.get('organism', {}).get('taxonId') == 9606 and key in {n.strip().casefold() for n in aliases}:
                    result['matches'].append({'accession': entry['primaryAccession'], 'names': aliases,
                        'organism': entry['organism'], 'url': 'https://www.uniprot.org/uniprotkb/' + entry['primaryAccession']})
            n = len(result['matches'])
            result['status'] = 'TRUNCATED' if len(entries) >= 100 else 'UNIQUE' if n == 1 else 'AMBIGUOUS' if n else 'NOT_FOUND_REVIEWED'
            self.failure_streak = 0
        except Exception as exc:
            result.update(status='REQUEST_FAILED', matches=[], error=type(exc).__name__ + ': ' + str(exc))
            self.failure_streak += 1
        self.entries[key] = result
        atomic_json(self.cache_path, {'schema': self.SCHEMA, 'entries': self.entries})
        return result

    def close(self):
        self.session.close()


UNIFIED_POLICY = 'unified-offline-cleaning-1.1.0'
# These are display fields only. Nested evidence, quotations, enums describing
# uncertainty, and negative observations (e.g. treatment="None") stay intact.
CLEAN_TEXT_FIELDS = {
    'gene_name', 'mutation_type', 'mutation_information', 'cell_line',
    'disease_information', 'literature_title', 'target_species', 'host_species',
    'normal_localization', 'mutated_localization', 'mislocalization_phenotype',
    'observed_protein',
}
CLEAN_MISSING = {'', 'unknown', 'not_reported', 'not reported', 'null', 'n/a', 'na', '*'}
CLEAN_FIELD_MAP = {
    'Gene Name': 'gene_name', 'Mutation Type': 'mutation_type',
    'Mutation Information': 'mutation_information',
    'Mislocalization Phenotype': 'mislocalization_phenotype',
    BASE_COLUMNS[6]: 'normal_localization', BASE_COLUMNS[7]: 'mutated_localization',
    'Cell Line': 'cell_line', 'Disease Information': 'disease_information',
    'Title': 'literature_title',
}
# Exact whole-label matches only. Mixed/chimeric/free-text labels stay pending.
# This routes extracted labels for review; it does not verify species from source.
CLEAN_NONHUMAN_LABELS = {
    'mouse', 'mus musculus', 'mouse (mus musculus)', 'murine',
    'rat', 'rattus norvegicus', 'chicken', 'gallus gallus', 'chicken (gallus gallus)',
    'arabidopsis thaliana', 'dictyostelium discoideum', 'xenopus', 'xenopus laevis',
    'xenopus tropicalis', 'medicago truncatula', 'escherichia coli',
    'stenotrophomonas maltophilia', 'pseudomonas aeruginosa', 'mycobacterium tuberculosis',
    'aspergillus nidulans', 'saccharomyces cerevisiae', 'drosophila melanogaster',
    'danio rerio', 'caenorhabditis elegans', 'bos taurus',
    'influenza a virus', 'influenza a virus a/pr/8/34', 'hepatitis b virus',
    'hepatitis c virus', 'hiv-1', 'human immunodeficiency virus 1',
    'human papillomavirus type 5', 'hpv5', 'herpes simplex virus 1', 'hsv-1',
}
SCREENING_LABELS = {
    'READY_FOR_SOURCE_REVIEW': '进入原文核验队列',
    'NEEDS_EVIDENCE': '待补证据',
    'OUT_OF_SCOPE_OR_DUPLICATE': '范围外或重复项',
}
SCREENING_REASONS = {
    'exact_duplicate': ('重复记录', '完整原始记录完全重复，保留并关联首条记录'),
    'explicit_nonhuman_target_label': ('物种范围', '抽取的目标物种标签明确非人源；范围初筛，仍保留原记录'),
    'target_species_not_explicitly_human': ('物种待核', '目标物种未报告、混合或有歧义，不根据宿主细胞推断'),
    'gene_identity_missing': ('身份缺失', '缺少明确基因身份'),
    'variant_identity_missing': ('身份缺失', '缺少明确变异身份'),
    'observed_entity_missing': ('测量对象', '缺少明确被观察蛋白名称'),
    'measured_entity_unresolved': ('测量对象', '突变蛋白与被观察对象的关系尚未明确或不同'),
    'localization_observation_missing': ('定位信息', '缺少可用的突变体定位或原文定位描述'),
    'outcome_unresolved': ('结果待核', '结果仍为未报告、不确定、冲突或未识别状态'),
    'quote_not_matched_to_source': ('引文匹配', '引文缺失、未匹配或缺少原文，需补充核对'),
    'reported_secondary_evidence': ('证据归属', '抽取标记为二手证据，原文核对后再决定是否排除'),
    'reported_downstream_effect': ('证据归属', '抽取标记为下游对象变化，原文核对后再决定是否排除'),
    'direct_localization_evidence_unconfirmed': ('实验方法', '直接定位证据类型尚未明确'),
    'extraction_conflict': ('信息冲突', '已有证据方法或结构化字段冲突，需要核对'),
    'no_change_comparator_missing': ('实验对照', '无变化结果缺少明确对照描述'),
}


def clean_display_value(value):
    """Map missing display values to null without interpreting 'None' as missing."""
    if value is None or isinstance(value, str) and value.strip().casefold() in CLEAN_MISSING:
        return None
    return value.strip() if isinstance(value, str) else value


def clean_record(record, rid, index, title, source_text, source_hash):
    """Pure, reversible formatting plus diagnostics; never a scientific verdict."""
    import copy
    if '_cleaning' in record:
        raise ValueError('Input is already cleaned. Re-run on the original extraction JSON.')
    result, changes = copy.deepcopy(record), []

    def change(key, value, rule):
        if result.get(key) != value:
            changes.append({'field': key, 'before': result.get(key), 'after': value, 'rule': rule,
                            'previously_present': key in result})
            result[key] = value

    change('pmid', rid.split(':')[0], 'numeric_pmid_as_text')
    for key in sorted(CLEAN_TEXT_FIELDS & record.keys()):
        if key == 'literature_title' and clean_display_value(title) is not None:
            continue
        change(key, clean_display_value(result[key]), 'trim_and_null_missing_display_value')
    if clean_display_value(title) is not None:
        change('literature_title', title.strip(), 'metadata_title_by_pmid')
    mutation, mutation_note = normalize_mutation(result.get('mutation_information'), result.get('mutation_type'))
    mutation_type, type_note = normalize_type(result.get('mutation_information'), result.get('mutation_type'))
    if result.get('mutation_information') is not None:
        change('mutation_information', clean_display_value(mutation), 'conservative_protein_notation')
    # No consequence is invented when both mutation fields are absent.
    if result.get('mutation_type') is not None or result.get('mutation_information') is not None:
        change('mutation_type', mutation_type, 'existing_mutation_type_rules')
    cell, cell_ids, cell_notes = normalize_cell(result.get('cell_line'))
    change('cell_line', clean_display_value(cell), 'unambiguous_cell_alias_only')
    if is_human_target(record):
        species_group = 'EXPLICIT_HUMAN'
        change('target_species', 'Homo sapiens', 'explicit_human_taxonomy_label')
    elif clean_display_value(record.get('target_species')) is None:
        species_group = 'SPECIES_UNREPORTED'
    else:
        species_group = 'OTHER_OR_AMBIGUOUS_SPECIES'

    quote = record.get('evidence_quote')
    if not source_text:
        quote_match = 'SOURCE_MISSING'
    elif clean_display_value(quote) is None:
        quote_match = 'QUOTE_MISSING'
    elif quote in source_text:
        quote_match = 'EXACT'
    elif re.sub(r'\s+', ' ', quote).strip() in re.sub(r'\s+', ' ', source_text).strip():
        quote_match = 'WHITESPACE_NORMALIZED'
    else:
        quote_match = 'NOT_FOUND'
    missing = [k for k in ['gene_name', 'mutation_information', 'normal_localization',
                           'mutated_localization'] if clean_display_value(result.get(k)) is None]
    scope_flags = []
    if (record.get('Evidence_Type') == 'SECONDARY_ONLY'
            or record.get('source_attribution_class') == 'SECONDARY_ONLY'
            or record.get('Outcome_Provenance') == 'SECONDARY_ONLY'):
        scope_flags.append('reported_secondary_evidence')
    if record.get('variant_observed_relationship') == 'DOWNSTREAM_EFFECT':
        scope_flags.append('reported_downstream_effect')
    if record.get('variant_observed_relationship') != 'SAME_PROTEIN':
        scope_flags.append('mutant_observed_identity_unresolved_or_different')
    row = {name: result.get(key) for name, key in CLEAN_FIELD_MAP.items()}
    row.update({name: clean_display_value(result.get(name)) for name in BASE_COLUMNS[16:]})
    flags = review_flags(row)
    conflicts = explicit_conflict_reasons(record, '; '.join(flags))
    result['_cleaning'] = {
        'policy': UNIFIED_POLICY, 'id': rid, 'original_record_index': index,
        'original_record_sha256': hashlib.sha256(json.dumps(record, ensure_ascii=False, sort_keys=True,
            separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest(),
        'source_sha256': source_hash, 'species_group': species_group,
        'quote_match': quote_match, 'missing_core_fields': missing,
        'reported_scope_flags': scope_flags, 'review_flags': flags, 'conflicts': conflicts,
        'mutation_note': mutation_note, 'type_note': type_note,
        'cell_ids': cell_ids, 'cell_notes': cell_notes, 'changes': changes,
        'duplicate_of': None, 'scientific_status': 'SOURCE_REVIEW_REQUIRED',
        'reference_sequence_verified': False, 'final_admission': False,
    }
    return result


def clean_export_row(record):
    """Keep the 26-column layout, retaining reported conditions as source wording."""
    row = {name: None for name in BASE_COLUMNS + [DETAIL]}
    row.update({name: record.get(key) for name, key in CLEAN_FIELD_MAP.items()})
    row.update({name: clean_display_value(record.get(name)) for name in BASE_COLUMNS[16:]})
    row.update(id=record['_cleaning']['id'], PMID=record['pmid'])
    details = []
    for key in ['outcome_class', 'normal_localization_raw', 'mutated_localization_raw',
                'mislocalization_phenotype_raw', 'observed_protein', 'comparator',
                'construct_background', 'treatment_condition', 'quantitative_result',
                'localization_event_type', 'evidence_location']:
        if clean_display_value(record.get(key)) is not None:
            details.append(key + ': ' + text(record[key]))
    row[DETAIL] = '\n'.join(details) or None
    return row


def screen_cleaned_record(record):
    """Assign one preliminary route. No records or experimental events are removed.

    Missing baseline/cell/tag/quantitation is not a blocking rule: a valid narrow
    observation can still proceed to source review. Extraction assertions about
    secondary/downstream evidence are flags, not automatic exclusion verdicts.
    """
    audit = record['_cleaning']; reasons, notes = [], []
    def usable(value):
        return (isinstance(value, str) and clean_display_value(value) is not None
                and value.strip().casefold() not in {'uncertain', 'unresolved', 'undetermined',
                    'other', 'none', 'not_applicable', 'not applicable', 'not available'})
    nonhuman = text(record.get('target_species')).strip().casefold() in CLEAN_NONHUMAN_LABELS
    if audit['duplicate_of']:
        reasons.append('exact_duplicate')
    if nonhuman:
        reasons.append('explicit_nonhuman_target_label')
    elif audit['species_group'] != 'EXPLICIT_HUMAN':
        reasons.append('target_species_not_explicitly_human')
    for field, reason in [('gene_name', 'gene_identity_missing'),
                          ('mutation_information', 'variant_identity_missing'),
                          ('observed_protein', 'observed_entity_missing')]:
        if not usable(record.get(field)):
            reasons.append(reason)
    if record.get('variant_observed_relationship') != 'SAME_PROTEIN':
        reasons.append('measured_entity_unresolved')
    if not any(usable(record.get(k)) for k in ['mutated_localization', 'mutated_localization_raw']):
        reasons.append('localization_observation_missing')
    if record.get('outcome_class') not in {'mislocalization', 'normal_experimentally_confirmed', 'observation_only'}:
        reasons.append('outcome_unresolved')
    if audit['quote_match'] not in {'EXACT', 'WHITESPACE_NORMALIZED'}:
        reasons.append('quote_not_matched_to_source')
    reasons.extend(flag for flag in audit['reported_scope_flags']
                   if flag in {'reported_secondary_evidence', 'reported_downstream_effect'})
    if record.get('Evidence_Type') not in {'SPATIAL_DIRECT', 'BIOCHEMICAL_LOCALIZATION'}:
        reasons.append('direct_localization_evidence_unconfirmed')
    if audit['conflicts']:
        reasons.append('extraction_conflict')
    if record.get('outcome_class') == 'normal_experimentally_confirmed' and not usable(record.get('comparator')):
        reasons.append('no_change_comparator_missing')
    if not usable(record.get('normal_localization')) and not usable(record.get('normal_localization_raw')):
        notes.append('缺少基线定位；可核对是否仅支持条件性观察，不直接排除')
    if not usable(record.get('comparator')):
        notes.append('缺少明确对照；不得直接宣称突变导致定位变化')
    if record.get('source_coverage') != 'COMPLETE':
        notes.append('原抽取未标记全文覆盖完整；核验时检查相关图表和实验条件')
    if not usable(record.get('cell_line')):
        notes.append('细胞名称缺失；可选信息不足不直接排除')
    status = ('OUT_OF_SCOPE_OR_DUPLICATE' if nonhuman or audit['duplicate_of'] else
              'NEEDS_EVIDENCE' if reasons else 'READY_FOR_SOURCE_REVIEW')
    audit.update(screening_status=status, screening_label=SCREENING_LABELS[status],
                 screening_reasons=reasons,
                 screening_reason_labels=[SCREENING_REASONS[r][1] for r in reasons],
                 problem_groups=list(dict.fromkeys(SCREENING_REASONS[r][0] for r in reasons)),
                 screening_notes=notes,
                 screening_basis='Extracted fields and deterministic diagnostics; not source adjudication')


def unified_clean(args):
    """One offline path from raw extraction to inspectable records and audit."""
    import csv
    from openpyxl import Workbook, load_workbook
    from io import BytesIO
    src, meta, md_dir, out = [p.resolve() for p in
        [args.input, args.metadata, args.md_dir, args.output_dir]]
    filenames = ['cleaned_records.json', 'cleaning_audit.xlsx', 'cleaning_report.json']
    if any(out / name in {src, meta} for name in filenames):
        raise ValueError('Output must not overwrite an input.')
    payload = json.loads(src.read_text(encoding='utf-8-sig'))
    records = payload.get('records') if isinstance(payload, dict) else payload
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError('Expected a JSON records array.')
    if any('_cleaning' in r for r in records):
        raise ValueError('Use the original extraction JSON, not a previous cleaning output.')
    with meta.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        if not {'PMID', 'Title'} <= set(reader.fieldnames or []):
            raise ValueError('Metadata requires PMID and Title columns.')
        metadata = {}
        for row in reader:
            pmid = row['PMID'].strip()
            if not re.fullmatch(r'[1-9]\d*', pmid) or pmid in metadata:
                raise ValueError('Invalid or duplicate metadata PMID: ' + pmid)
            metadata[pmid] = row
    sources, source_hashes = {}, {}
    for pmid in sorted(metadata):
        path = md_dir / (pmid + '.md')
        content = path.read_bytes() if path.exists() else None
        sources[pmid] = content.decode('utf-8-sig') if content is not None else ''
        source_hashes[pmid] = hashlib.sha256(content).hexdigest() if content is not None else None
    binding = {'policy': UNIFIED_POLICY, 'input_sha256': file_sha256(src),
               'metadata_sha256': file_sha256(meta), 'source_sha256': source_hashes,
               'script_sha256': file_sha256(Path(__file__))}
    previous = out / 'cleaning_report.json'
    if any((out / name).exists() for name in filenames) and not args.force:
        if not previous.exists() or json.loads(previous.read_text(encoding='utf-8'))['binding'] != binding:
            raise ValueError('Inputs or rules changed; choose a new output directory or use --force.')
    cleaned, counts, exact_seen = [], Counter(), {}
    for index, record in enumerate(records, 1):
        pmid = text(record.get('pmid')).strip()
        if pmid not in metadata:
            raise ValueError('Record PMID missing from metadata: ' + pmid)
        counts[pmid] += 1
        item = clean_record(record, pmid + ':' + str(counts[pmid]).zfill(4), index,
                            metadata[pmid]['Title'], sources[pmid], source_hashes[pmid])
        audit = item['_cleaning']; fingerprint = audit['original_record_sha256']
        audit['duplicate_of'] = exact_seen.get(fingerprint)
        exact_seen.setdefault(fingerprint, audit['id'])
        screen_cleaned_record(item)
        cleaned.append(item)
    report = {
        'policy': UNIFIED_POLICY, 'binding': binding, 'input_records': len(records),
        'output_records': len(cleaned), 'dropped_records': 0, 'papers': len(metadata),
        'papers_with_records': len(counts), 'papers_without_records': len(metadata) - len(counts),
        'species_groups': dict(Counter(r['_cleaning']['species_group'] for r in cleaned)),
        'quote_matches': dict(Counter(r['_cleaning']['quote_match'] for r in cleaned)),
        'changed_records': sum(bool(r['_cleaning']['changes']) for r in cleaned),
        'field_changes': dict(Counter(c['field'] for r in cleaned for c in r['_cleaning']['changes'])),
        'exact_duplicates_retained': sum(r['_cleaning']['duplicate_of'] is not None for r in cleaned),
        'records_with_missing_core_fields': sum(bool(r['_cleaning']['missing_core_fields']) for r in cleaned),
        'records_with_conflicts': sum(bool(r['_cleaning']['conflicts']) for r in cleaned),
        'screening_counts': {status: sum(r['_cleaning']['screening_status'] == status for r in cleaned)
                             for status in SCREENING_LABELS},
        'human_screening_counts': {status: sum(r['_cleaning']['species_group'] == 'EXPLICIT_HUMAN'
                                     and r['_cleaning']['screening_status'] == status for r in cleaned)
                                   for status in SCREENING_LABELS},
        'problem_group_counts': dict(Counter(g for r in cleaned for g in r['_cleaning']['problem_groups'])),
        'screening_reason_counts': dict(Counter(g for r in cleaned for g in r['_cleaning']['screening_reasons'])),
        'scientific_status': 'SOURCE_REVIEW_REQUIRED', 'new_api_requests': 0,
        'network_enabled': False, 'final_admission': False, 'columns': 26,
    }
    columns = list(BASE_COLUMNS); columns.insert(columns.index('Mislocalization Phenotype') + 1, DETAIL)
    rows = [clean_export_row(r) for r in cleaned]
    book = Workbook(write_only=True); expected = {}
    def sheet(name, headers, values, widths=None):
        expected[name] = write_sheet(book, name, headers, values, widths)
    sheet('Summary', ['Metric', 'Value'], [[k, v] for k, v in report.items() if k != 'binding'] + [
        ['interpretation', 'All input records retained. Cleaning and quote matching do not establish scientific validity.'],
        ['species_groups', 'Reported target species only; host cells never establish protein species.'],
        ['missingness', 'Display missing markers become null. None/no treatment and uncertainty enums are retained.'],
        ['duplicates', 'Only identical complete raw records are linked; no experimental events are merged.'],
        ['evidence', 'Quotations and nested evidence are unchanged; normalized values require renewed source review.'],
        ['screening', 'Three mutually exclusive preliminary routes; READY means ready for review, not accepted or A/B.'],
        ['scope_basis', 'Exact nonhuman target labels and complete raw duplicates are routed aside but retained. Ambiguous species, secondary and downstream flags require evidence review.'],
        ['optional_fields', 'Missing baseline, cell line, tag or quantitation does not by itself block entry to source review.'],
        ['problem_groups', 'Problem counts overlap; only screening_counts are mutually exclusive.'],
        ['replay', 'Run again on the original extraction. Matching inputs and rules may reuse this output directory.'],
    ], [43, 110])
    sheet('Records', columns, [[row[h] for h in columns] for row in rows])
    sheet('Audit', ['id', 'PMID', '初筛去向', '初筛原因', '问题分组', '核验提示', 'Screening code',
                   'Species group', 'Reported target species', 'Host species',
                   'Quote match', 'Missing core fields', 'Reported scope flags', 'Conflicts',
                   'Review flags', 'Exact duplicate of', 'Scientific status', 'Mutation note',
                   'Type note', 'Cell notes', 'Source SHA256', 'Original record SHA256',
                   'Input error flags', 'Reported outcome', 'Reported evidence type',
                   'Reported relationship', 'Reported source coverage'], [
        [r['_cleaning']['id'], r['pmid'], r['_cleaning']['screening_label'],
         '\n'.join(r['_cleaning']['screening_reason_labels']) or '基本信息可用，进入原文核验；尚未收录',
         '; '.join(r['_cleaning']['problem_groups']), '\n'.join(r['_cleaning']['screening_notes']),
         r['_cleaning']['screening_status'], r['_cleaning']['species_group'], original.get('target_species'),
         original.get('host_species')] + [r['_cleaning'][k] for k in ['quote_match', 'missing_core_fields',
         'reported_scope_flags', 'conflicts', 'review_flags', 'duplicate_of', 'scientific_status',
         'mutation_note', 'type_note', 'cell_notes', 'source_sha256', 'original_record_sha256']] +
        [original.get(k) for k in ['error_type', 'outcome_class', 'Evidence_Type',
                                   'variant_observed_relationship', 'source_coverage']]
        for r, original in zip(cleaned, records)])
    sheet('Changes', ['id', 'PMID', 'Field', 'Original', 'Cleaned', 'Rule'], [
        [r['_cleaning']['id'], r['pmid'], c['field'], c['before'], c['after'], c['rule']]
        for r in cleaned for c in r['_cleaning']['changes']])
    sheet('Papers', ['PMID', 'Title', 'Records', 'Source available', 'Source SHA256'], [
        [pmid, p['Title'], counts[pmid], bool(sources[pmid]), source_hashes[pmid]]
        for pmid, p in metadata.items()], [16, 100, 16, 20, 70])
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='cleaning_', dir=out) as staging:
        staging = Path(staging)
        workbook = staging / 'cleaning_audit.xlsx'; book.save(workbook)
        check = load_workbook(BytesIO(workbook.read_bytes()), read_only=True, data_only=False)
        try:
            for name, digest in expected.items():
                if data_digest(check[name].values) != digest:
                    raise ValueError('Workbook readback mismatch: ' + name)
        finally:
            check.close()
        for path, digest in [(src, binding['input_sha256']), (meta, binding['metadata_sha256'])]:
            if file_sha256(path) != digest:
                raise ValueError('Input changed during cleaning.')
        for pmid, digest in source_hashes.items():
            path = md_dir / (pmid + '.md')
            if (file_sha256(path) if path.exists() else None) != digest:
                raise ValueError('Source changed during cleaning: ' + pmid)
        report['records_data_sha256'] = data_digest([[row[h] for h in columns] for row in rows])
        report['readback_verified'] = True
        atomic_json(staging / 'cleaned_records.json', {'schema_version': UNIFIED_POLICY,
                    'input_sha256': binding['input_sha256'], 'records': cleaned})
        report['cleaned_json_sha256'] = file_sha256(staging / 'cleaned_records.json')
        report['workbook_sha256'] = file_sha256(workbook)
        atomic_json(staging / 'cleaning_report.json', report)
        for name in filenames:
            os.replace(staging / name, out / name)
    print(json.dumps({k: v for k, v in report.items() if k != 'binding'}, ensure_ascii=False, indent=2))
    return report


def unified_clean_main():
    cli = argparse.ArgumentParser(description='One offline cleaning step; retain every record and audit every change.')
    cli.add_argument('--unified', action='store_true')
    cli.add_argument('--preset', choices=['original', 'expanded'], default='original')
    cli.add_argument('--input', type=Path)
    cli.add_argument('--metadata', type=Path)
    cli.add_argument('--md-dir', type=Path)
    cli.add_argument('--output-dir', type=Path)
    cli.add_argument('--force', action='store_true', help='Replace this cleaning output after inputs or rules change')
    cli.add_argument('--version', action='version', version=UNIFIED_POLICY)
    args = cli.parse_args()
    args.input = args.input or ROOT / 'data/04_extract_deepseek' / args.preset / 'extracted_records.json'
    args.metadata = args.metadata or ROOT / 'data/01_get_pubmed_metadata' / ('papers_' + args.preset + '.csv')
    args.md_dir = args.md_dir or ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean'
    args.output_dir = args.output_dir or ROOT / 'data/05_clean_records' / args.preset
    unified_clean(args)
    return 0


def summarize_main():
    if '--unified' in sys.argv[1:]:
        return unified_clean_main()
    import csv
    from openpyxl import Workbook, load_workbook
    cli = argparse.ArgumentParser(description='Summarize Step 04 JSON using legacy-table normalization rules 1.1.0.')
    cli.add_argument('--preset', choices=['original', 'expanded'], default='original')
    cli.add_argument('--input', type=Path)
    cli.add_argument('--metadata', type=Path, help='Paper metadata CSV with PMID and Title columns')
    cli.add_argument('--output-dir', type=Path)
    cli.add_argument('--offline', action='store_true', help='Use metadata-table titles without network requests (default)')
    cli.add_argument('--online', action='store_true', help='Verify titles online using the 1.1.0 resolver')
    cli.add_argument('--cache', type=Path, help='Optional title-verification cache used only with --online')
    cli.add_argument('--uniprot', action='store_true', help='Optionally verify human base-protein names via UniProt; never alter localization or reserved fields')
    cli.add_argument('--uniprot-cache', type=Path, help='UniProt cache path; default: output directory/uniprot_cache.json')
    cli.add_argument('--md-dir', type=Path, help='Source Markdown for evidence checks')
    cli.add_argument('--quality', choices=['strict', 'basic'], default='strict', help='Default strict: require source-supported identity and core evidence')
    cli.add_argument('--force', action='store_true')
    cli.add_argument('--complete-only', action='store_true', help='Keep records with complete key fields; retain exclusions and reasons in separate sheets')
    cli.add_argument('--self-test', action='store_true')
    cli.add_argument('--version', action='version', version='1.1.0-json-adapter.2')
    args = cli.parse_args()
    if args.self_test:
        return 0 if self_test() else 2
    if args.online and args.offline:
        cli.error('--online and --offline are mutually exclusive')
    if args.uniprot and args.offline:
        cli.error('--uniprot requires network access; omit --offline. Titles remain offline without --online.')
    src = (args.input or ROOT / 'data/04_extract_deepseek' / args.preset / 'extracted_records.json').resolve()
    meta = (args.metadata or ROOT / 'data/01_get_pubmed_metadata' / ('papers_' + args.preset + '.csv')).resolve()
    out = (args.output_dir or ROOT / 'data/05_summarize_records' / args.preset).resolve()
    md_dir = (args.md_dir or ROOT / 'data/03_mineru_pdf_to_md' / args.preset / 'clean').resolve()
    target = out / ('summary_clean.xlsx' if args.complete_only else 'summary.xlsx')
    lookup_cache = (args.uniprot_cache or out / 'uniprot_cache.json').resolve()
    for cache in [args.cache, lookup_cache if args.uniprot else None]:
        if cache and cache.resolve() in {src, meta, target}:
            cli.error('Cache must not replace inputs or workbook output.')
    if target.exists() and not args.force:
        cli.error('Output exists; use another output directory or --force.')
    input_hash, metadata_hash = file_sha256(src), file_sha256(meta)
    payload = json.loads(src.read_text(encoding='utf-8-sig'))
    records = payload.get('records') if isinstance(payload, dict) else payload
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError('Expected a JSON records array.')
    with meta.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        if not {'PMID', 'Title'} <= set(reader.fieldnames or []):
            raise ValueError('Metadata requires PMID and Title headers.')
        papers = list(reader)
    metadata = {}
    for paper in papers:
        pmid = text(paper['PMID']).strip()
        if not re.fullmatch(r'[1-9]\d*', pmid) or pmid in metadata:
            raise ValueError('Invalid or duplicate metadata PMID: ' + pmid)
        metadata[pmid] = paper
    field_map = {'Gene Name': 'gene_name', 'Mutation Type': 'mutation_type',
                 'Mutation Information': 'mutation_information',
                 'Mislocalization Phenotype': 'mislocalization_phenotype',
                 BASE_COLUMNS[6]: 'normal_localization', BASE_COLUMNS[7]: 'mutated_localization',
                 'Cell Line': 'cell_line', 'Disease Information': 'disease_information'}
    source_rows, auxiliary, provenance, raw_chunks = [], {}, [], []
    per_paper = Counter()
    for index, record in enumerate(records, 1):
        pmid = text(record.get('pmid')).strip()
        if pmid not in metadata:
            raise ValueError('Record PMID missing from metadata: ' + pmid)
        per_paper[pmid] += 1
        rid = pmid + ':' + str(per_paper[pmid]).zfill(4)
        row = {h: None for h in BASE_COLUMNS}
        row.update({name: record.get(key) for name, key in field_map.items()})
        row.update({name: record.get(name) for name in BASE_COLUMNS[16:]})
        row.update({'id': rid, 'PMID': pmid,
                    'Title': metadata[pmid].get('Title') or record.get('literature_title')})
        source_rows.append((index + 1, row))
        auxiliary[(rid, pmid)] = record
        provenance.append([rid, pmid, index, 'records[' + str(index - 1) + ']',
                           record.get('candidate_id'), record.get('literature_title'), row['Title'],
                           'Metadata CSV by PMID; external verification depends on title mode'])
        serialized = json.dumps(record, ensure_ascii=False, separators=(',', ':'))
        for offset in range(0, len(serialized), 28000):
            raw_chunks.append([rid, pmid, offset // 28000 + 1, serialized[offset:offset + 28000]])
    # Metadata titles are authoritative for the default offline workflow.
    # Cached or remote titles can replace them only when --online is explicit.
    resolver = TitleResolver('online' if args.online else 'keep', args.cache if args.online else None)
    resolver.resolve(set(metadata))
    headers, rows, changes, statuses, mappings, counts = normalize_rows(
        BASE_COLUMNS, source_rows, auxiliary, set(), resolver)
    assert len(rows) == len(records) and len(headers) == 26
    assert all(all(r[h] is None for h in BLANK_COLUMNS) for r in rows)
    counts['title_unverified_pmids'] = sum(not resolver.usable(p) for p in metadata)
    counts['aux_loaded_keys'] = len(auxiliary)
    report = {'adapter_version': '1.1.0-json-adapter.2', 'normalizer_version': VERSION,
              'normalizer_source_sha256': CORE_SOURCE_SHA256,
              'input_sha256': input_hash, 'metadata_sha256': metadata_hash,
              'script_sha256': file_sha256(Path(__file__)),
              'papers': len(metadata), 'papers_with_records': len(per_paper),
              'papers_without_records': len(metadata) - len(per_paper),
              'records': len(records), 'columns': len(headers), 'title_mode': resolver.mode,
              'counts': dict(counts)}
    required = ['PMID', 'Gene Name', 'Mutation Type', 'Mutation Information',
                'Mislocalization Phenotype', BASE_COLUMNS[6], BASE_COLUMNS[7],
                'Cell Line', 'Title', 'Experimental_Method', 'Evidence_Type', 'Outcome_Provenance']
    missing_markers = {'', 'unknown', 'not_reported', 'undetermined', 'null', 'n/a', 'na', '*', 'unresolved'}
    def missing_key(value):
        value = text(value).strip().casefold()
        return value in missing_markers or any(part.strip() in missing_markers
            for part in re.split(r'[>;]', value))
    retained, excluded, species_rows = [], [], []
    normalization_flags = {(s[1], s[2]): s[-1] for s in statuses}
    conflict_counts = Counter()
    conflict_record_count = 0
    for row in rows:
        original = auxiliary[(row['id'], row['PMID'])]
        human = is_human_target(original)
        species_rows.append([row['id'], row['PMID'], original.get('target_species'),
                             original.get('host_species'), human])
        absent = [name for name in required if missing_key(row.get(name))] if args.complete_only else []
        if args.complete_only:
            conflicts = explicit_conflict_reasons(original, normalization_flags.get((row['id'], row['PMID']), ''))
            absent.extend(conflicts)
            conflict_counts.update(conflicts)
            conflict_record_count += bool(conflicts)
        if not human:
            absent.append('target_species_not_explicitly_human')
        if absent:
            excluded.append([row[h] for h in headers] + ['; '.join(absent)])
        else:
            retained.append(row)
    shown = retained
    report.update({'completeness_filter': args.complete_only, 'retained_records': len(shown),
                   'excluded_records': len(excluded), 'human_only': True,
                   'explicit_human_input_records': sum(r[-1] for r in species_rows),
                   'retained_papers': len({r['PMID'] for r in shown}),
                   'conflict_records_excluded': conflict_record_count,
                   'conflict_reasons': dict(conflict_counts),
                   'required_fields': required if args.complete_only else []})
    lookup_rows = []
    lookup_counts = Counter()
    if args.uniprot:
        lookup = HumanUniProtLookup(lookup_cache)
        try:
            for i, row in enumerate(shown, 1):
                original = auxiliary[(row['id'], row['PMID'])]
                a = lookup.lookup(original.get('gene_name'))
                b = lookup.lookup(original.get('observed_protein'))
                left = a['matches'][0]['accession'] if a['status']=='UNIQUE' else ''
                right = b['matches'][0]['accession'] if b['status']=='UNIQUE' else ''
                mapping = ('SAME_BASE_PROTEIN' if left==right else 'DISTINCT_BASE_PROTEINS') if left and right else 'SKIPPED'
                lookup_counts[mapping] += 1
                lookup_rows.append([row['id'], row['PMID'], original.get('gene_name'),
                    original.get('observed_protein'), original.get('target_species'),
                    a['status'], b['status'], left, right, mapping,
                    a.get('queried_at'), b.get('queried_at'),
                    a.get('response_url'), b.get('response_url'),
                    original.get('variant_observed_relationship'), original.get('relationship_verification')])
                if i == 1 or i % 10 == 0 or i == len(shown):
                    print(f'UniProt: {i}/{len(shown)} records; {lookup.requests} requests, {lookup.cache_hits} cache hits', flush=True)
            report.update(uniprot_name_mapping=dict(lookup_counts), uniprot_requests=lookup.requests,
                          uniprot_cache_hits=lookup.cache_hits, uniprot_cache_sha256=file_sha256(lookup_cache) if lookup_cache.exists() else None)
        finally:
            lookup.close()
    quality_rows, pending_rows = [], []
    basic_count = len(shown)
    if args.quality == 'strict':
        high_quality = []
        lookup_mapping = {r[0]:r[9] for r in lookup_rows}
        sources, source_hashes = {}, {}
        for row in shown:
            pmid = row['PMID']
            if pmid not in sources:
                md_file = md_dir / (pmid + '.md')
                sources[pmid] = md_file.read_text(encoding='utf-8-sig') if md_file.exists() else ''
                source_hashes[pmid] = file_sha256(md_file) if md_file.exists() else None
            original = auxiliary[(row['id'], pmid)]
            category, reasons = strict_quality_decision(original, sources[pmid], lookup_mapping.get(row['id'], 'SKIPPED'))
            quality_rows.append([row['id'], pmid, category, '; '.join(reasons)])
            values = [row[h] for h in headers]
            if category == 'HIGH_CONFIDENCE_CANDIDATE':
                high_quality.append(row)
            elif category == 'NEEDS_REVIEW':
                pending_rows.append(values + ['; '.join(reasons)])
            else:
                excluded.append(values + ['; '.join(reasons)])
        shown = high_quality
        report['source_markdown_sha256'] = source_hashes
    report.update(quality_mode=args.quality, basic_candidate_records=basic_count,
                  retained_records=len(shown), review_records=len(pending_rows), excluded_records=len(excluded),
                  retained_papers=len({r['PMID'] for r in shown}),
                  review_reason_counts=dict(Counter(reason for r in pending_rows for reason in r[-1].split('; '))))
    review_groups = defaultdict(list)
    for values in pending_rows:
        review_groups[values[headers.index('PMID')]].append(values)
    report['review_papers'] = len(review_groups)
    assert len(shown) + len(pending_rows) + len(excluded) == len(records)
    summary = [[k, v] for k, v in report.items() if k not in {'counts', 'source_markdown_sha256'}]
    summary += [[k, v] for k, v in sorted(counts.items())]
    summary += [['scope', 'Explicit human target species only; exclusions and all original records preserved'],
                ['uniprot_scope', 'Reviewed human entries, exact names, unique matches only; base-protein identity does not confirm mutant identity'],
                ['title_source', 'Metadata CSV by PMID; original extracted title used if metadata title is empty. External verification requires --online.'],
                ['raw_json', 'Concatenate Raw JSON chunks by id and Chunk_Order to reconstruct each record'],
                ['interpretation', 'Formatting normalization is not scientific validation']]
    out.mkdir(parents=True, exist_ok=True)
    if args.cache and args.cache.resolve() in {src, meta, target}:
        raise ValueError('Cache must not replace inputs or output.')
    book = Workbook(write_only=True)
    expected = {}
    def sheet(name, fields, values, widths=None):
        expected[name] = write_sheet(book, name, fields, values, widths)
    sheet('Summary', ['Metric', 'Value'], summary, [55, 110])
    sheet('Normalized Records', headers, ([r[h] for h in headers] for r in shown))
    sheet('Excluded Records', headers + ['Exclusion_Reasons'], excluded)
    if args.quality == 'strict':
        sheet('Needs Review', headers + ['Review_Reasons'], pending_rows)
        sheet('Quality Audit', ['id', 'PMID', 'Category', 'Reasons'], quality_rows)
        sheet('Review by Paper', ['PMID', 'Title', 'Pending_Records', 'Record_IDs', 'Reason_Counts'],
              [[pmid, metadata[pmid].get('Title'), len(values), '; '.join(r[0] for r in values),
                dict(Counter(reason for r in values for reason in r[-1].split('; ')))]
               for pmid, values in review_groups.items()], [16, 85, 18, 60, 100])
    sheet('Species Audit', ['id', 'PMID', 'Target_Species', 'Host_Species', 'Explicit_Human_Target'], species_rows)
    if args.uniprot:
        sheet('Protein Lookup', ['id', 'PMID', 'Gene_Name', 'Observed_Protein', 'Target_Species',
              'Gene_Lookup_Status', 'Observed_Lookup_Status', 'Gene_UniProt', 'Observed_UniProt',
              'Base_Protein_Mapping', 'Gene_Query_UTC', 'Observed_Query_UTC', 'Gene_Query_URL',
              'Observed_Query_URL', 'Extracted_Relationship', 'Relationship_Verification'], lookup_rows)
    retained_per_paper = Counter(r['PMID'] for r in shown)
    sheet('Papers', ['PMID', 'Title', 'Input_Records', 'Retained_Records'],
          [[p, r.get('Title'), per_paper[p], retained_per_paper[p]] for p, r in metadata.items()], [16, 100, 18, 18])
    sheet('Original Records', BASE_COLUMNS, ([r[h] for h in BASE_COLUMNS] for _, r in source_rows))
    sheet('Field Changes', ['Source_Row', 'id', 'PMID', 'Field', 'Original', 'Normalized', 'Reason'], changes)
    sheet('Review Status', ['Source_Row', 'id', 'PMID', 'HGVS_Status', 'Localization_Basis',
                           'Cellosaurus_ID', 'Title_Status', 'Review_Flags'], statuses)
    sheet('Type Mapping', ['Original_Type', 'Normalized_Type', 'Records', 'Reason'],
          [[a, b, n, why] for (a, b, why), n in sorted(mappings.items())])
    sheet('Provenance', ['id', 'PMID', 'Record_Index', 'JSON_Locator', 'Candidate_ID',
                         'Extracted_Title', 'Metadata_Title', 'Title_Basis'], provenance)
    sheet('Raw JSON', ['id', 'PMID', 'Chunk_Order', 'JSON_Text'], raw_chunks, [22, 16, 16, 110])
    sheet('Rules', ['Item', 'Rule', 'Reference'], make_rules(counts, resolver, 'JSON records', []))
    handle, temporary = tempfile.mkstemp(suffix='.xlsx', dir=out)
    os.close(handle)
    try:
        book.save(temporary)
        check = load_workbook(temporary, read_only=True, data_only=False)
        try:
            for name, digest in expected.items():
                if data_digest(check[name].values) != digest:
                    raise ValueError('Readback mismatch: ' + name)
        finally:
            check.close()
        if file_sha256(src) != input_hash or file_sha256(meta) != metadata_hash:
            raise ValueError('Input changed during execution.')
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)
    atomic_json(out / ('summary_clean_report.json' if args.complete_only else 'summary_report.json'), report)
    print(json.dumps({'output': str(target), 'papers': len(metadata), 'records': len(records),
                      'retained': len(shown), 'needs_review': len(pending_rows), 'excluded': len(excluded),
                      'columns': len(headers), 'title_mode': resolver.mode}, indent=2))
    return 0


