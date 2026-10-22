#!/usr/bin/env python3
"""Fetch HIV-1 reference genomes + a diverse strain set from NCBI (E-utilities).

Run with /usr/bin/python3 (has biopython). Outputs:
  data/hiv/M19921_NL4-3.fasta   (provisional: proviral slice extracted here in-place)
  data/hiv/strains/<acc>.fasta  (diverse complete genomes)
  data/hiv/strains.tsv          (accession, title, length, subtype label, country)

Idempotent: re-running skips already-fetched items. ~3 req/s to NCBI (no key).
"""
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HIV_DIR = os.path.join(HERE, "data", "hiv")
STRAINS_DIR = os.path.join(HIV_DIR, "strains")
TOOL = "crispr-hiv-ward"


def _get(url, tries=6):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": TOOL})
            with urllib.request.urlopen(req, timeout=60) as r:
                out = r.read().decode("utf-8", "replace")
                time.sleep(0.5)  # global pacing: stay under 3 req/s
                return out
        except urllib.error.HTTPError as e:
            if e.code == 429 and i < tries - 1:
                time.sleep(min(20, 2 ** (i + 1)))
                continue
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                raise
            print(f"  retry {i+1}: {e}", file=sys.stderr)
            time.sleep(2 * (i + 1))


def esearch(term, retmax=500):
    q = urllib.parse.urlencode({"db": "nuccore", "retmode": "json", "retmax": retmax,
                                "term": term, "tool": TOOL})
    d = json.loads(_get(EUTILS + "esearch.fcgi?" + q))
    return d["esearchresult"]["idlist"]


def esummary(ids):
    out = {}
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        q = urllib.parse.urlencode({"db": "nuccore", "retmode": "json",
                                    "id": ",".join(chunk), "tool": TOOL})
        d = json.loads(_get(EUTILS + "esummary.fcgi?" + q))
        res = d.get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
        time.sleep(0.4)
    return out


def efetch_fasta(ids):
    parts = []
    for i in range(0, len(ids), 20):
        q = urllib.parse.urlencode({"db": "nuccore", "rettype": "fasta",
                                    "retmode": "text", "id": ",".join(ids[i:i + 20]),
                                    "tool": TOOL})
        parts.append(_get(EUTILS + "efetch.fcgi?" + q))
        time.sleep(0.4)
    return "".join(parts)


def parse_fasta(text):
    """Yield (name, seq) from FASTA text."""
    name, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def parse_subtype_country(title):
    sub = None
    m = re.search(r"(?:subtype|clade)\s*[-_ ]?([A-K]|CRF\d+_?[A-Z]{0,2})", title, re.I)
    if m:
        sub = m.group(1).upper().replace("_", "_")
    m2 = re.search(r"\bfrom\s+([A-Za-z ]+?)(?:,|;|\(|\b\d{4}\b)", title)
    country = m2.group(1).strip() if m2 else ""
    return sub, country


def extract_nl43():
    """M19921 is the pNL4-3 plasmid (14.8 kb with bacterial backbone).
    Slice the HIV-1 proviral span from its GenBank annotation."""
    from Bio import SeqIO
    gb_path = os.path.join(HIV_DIR, "M19921_NL4-3.gb")
    if not os.path.exists(gb_path):
        q = urllib.parse.urlencode({"db": "nuccore", "rettype": "gb",
                                    "retmode": "text", "id": "M19921", "tool": TOOL})
        open(gb_path, "w").write(_get(EUTILS + "efetch.fcgi?" + q))
        time.sleep(0.4)
    rec = SeqIO.read(gb_path, "genbank")
    spans = []
    for f in rec.features:
        note = f.qualifiers.get("note", [""])[0]
        if f.type == "repeat_region" and "LTR" in note:
            spans.append((int(f.location.start), int(f.location.end)))
    if len(spans) < 2:
        raise RuntimeError("Expected 5' and 3' LTR features in M19921")
    lo, hi = min(s for s, _ in spans), max(e for _, e in spans)
    prov = rec.seq[lo:hi]
    print(f"NL4-3 proviral span: {lo}-{hi} -> {len(prov)} bp")
    if not 8500 <= len(prov) <= 9900:
        raise RuntimeError(f"Unexpected proviral length {len(prov)}")
    out = os.path.join(HIV_DIR, "NL4-3_proval.fasta")
    with open(out, "w") as fh:
        fh.write(">NL4-3_provirus\n" + str(prov).upper() + "\n")
    return out


def fetch_strain_set(target=24, per_subtype=3):
    os.makedirs(STRAINS_DIR, exist_ok=True)
    term = ('Human immunodeficiency virus type 1[Organism] AND '
            'complete genome[Title] AND 8500[SLEN]:9900[SLEN]')
    ids = esearch(term, 2000)
    print(f"esearch found {len(ids)} complete HIV-1 genomes")
    summ = esummary(ids)
    labeled, unlabeled = [], []
    for uid, s in summ.items():
        slen = s.get("slen", 0)
        acc = s.get("accessionversion", uid)
        if not 8500 <= slen <= 9900:
            continue
        title = s.get("title", "")
        if any(w in title.lower() for w in ("mutant", "synthetic", "vector")):
            continue  # engineered constructs would bias conservation
        sub, country = parse_subtype_country(title)
        r = {"uid": uid, "acc": acc, "title": title,
             "len": slen, "subtype": sub or "", "country": country}
        (labeled if sub else unlabeled).append(r)
    picked, seen_sub, seen_country = [], {}, {}

    def push(r):
        if len(picked) >= target:
            return False
        c = (r["country"] or "?").lower()
        if r["subtype"]:
            cap = 5 if r["subtype"] == "B" else per_subtype
            if seen_sub.get(r["subtype"], 0) >= cap:
                return False
            seen_sub[r["subtype"]] = seen_sub.get(r["subtype"], 0) + 1
        if seen_country.get(c, 0) >= 3:
            return False
        seen_country[c] = seen_country.get(c, 0) + 1
        picked.append(r)
        return True

    for r in sorted(labeled, key=lambda x: (x["subtype"] != "B", x["subtype"])):
        push(r)
    for r in sorted(unlabeled, key=lambda x: -x["len"]):
        push(r)
    print(f"selected {len(picked)} strains; subtype labels: {seen_sub}")
    if len(picked) < 12:
        raise RuntimeError(f"too few strains selected ({len(picked)})")
    fa = efetch_fasta([r["uid"] for r in picked])
    by_name = {name: seq for name, seq in parse_fasta(fa)}
    rows = []
    for r in picked:
        # match fetched record by accession prefix
        seq = None
        for name, s in by_name.items():
            if name.startswith(r["acc"].split(".")[0]):
                seq = s
                break
        if seq is None or not 8500 <= len(seq) <= 9900:
            print(f"  SKIP {r['acc']} (no/odd sequence)")
            continue
        path = os.path.join(STRAINS_DIR, r["acc"] + ".fasta")
        with open(path, "w") as fh:
            fh.write(f">{r['acc']}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
        rows.append((r["acc"], r["title"][:90], len(seq),
                     r["subtype"] or "", r["country"]))
    with open(os.path.join(HIV_DIR, "strains.tsv"), "w") as fh:
        fh.write("accession\ttitle\tlength\tsubtype\tcountry\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"wrote {len(rows)} strain files + strains.tsv")


def main():
    os.makedirs(HIV_DIR, exist_ok=True)
    if not os.path.exists(os.path.join(HIV_DIR, "NL4-3_proval.fasta")):
        print("Extracting NL4-3 proviral genome from M19921...")
        extract_nl43()
    strains_tsv = os.path.join(HIV_DIR, "strains.tsv")
    if not os.path.exists(strains_tsv):
        fetch_strain_set()
    with open(strains_tsv) as fh:
        n = sum(1 for line in fh if not line.startswith("accession"))
    print(f"strain set ready: {n} strains")


if __name__ == "__main__":
    main()