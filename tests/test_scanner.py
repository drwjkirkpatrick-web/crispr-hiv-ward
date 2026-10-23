#!/usr/bin/env python3
"""crispr_scan test suite — plain python3 asserts (no pytest).

Run:  python3 tests/test_scanner.py     (or `make -C src test`)
Builds ../src/crispr_scan via make if missing.  Scratch files are written
under TMPDIR and cleaned up; nothing is written inside the repo.
"""

import contextlib
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(ROOT, "..", "src")
BIN = os.path.join(SRC_DIR, "crispr_scan")
GRCH38_BP = 3_100_000_000

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------


def ensure_bin():
    if not os.path.exists(BIN):
        r = subprocess.run(["make", "-C", SRC_DIR], capture_output=True, text=True)
        assert r.returncode == 0, f"make failed:\n{r.stdout}\n{r.stderr}"


def run(*args, expect_ok=True, timeout=600):
    cmd = [BIN] + [str(a) for a in args]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if expect_ok:
        assert r.returncode == 0, (
            f"cmd {cmd} rc={r.returncode}\nstdout={r.stdout!r}\nstderr={r.stderr!r}"
        )
    return r


@contextlib.contextmanager
def tmpdir():
    with tempfile.TemporaryDirectory(prefix="crispr_scan_test_") as d:
        yield d


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, idxs):
    """flip the bases at 0-based idxs to a different base"""
    out = list(s)
    for i in idxs:
        out[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[out[i]]
    return "".join(out)


def write_fasta(path, records, eol="\n", width=None, blank_lines=False):
    """records: list of (header_line, seq).  Multi-line if width given."""
    parts = []
    for hdr, seq in records:
        parts.append(">" + hdr)
        if width:
            for i in range(0, len(seq), width):
                parts.append(seq[i : i + width])
        else:
            parts.append(seq)
        if blank_lines:
            parts.append("")
    with open(path, "wb") as f:
        f.write(eol.join(parts).encode() + eol.encode())


def read_fasta(path):
    recs = []
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                toks = line[1:].split()
                name = toks[0] if toks else ""
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_guides(path, guides, header=True):
    with open(path, "w") as f:
        if header:
            f.write("#guide_id\tsystem\tsequence\n")
        for gid, sysname, seq in guides:
            f.write(f"{gid}\t{sysname}\t{seq}\n")


def parse_tsv(path):
    with open(path) as f:
        lines = f.read().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln]
    return header, rows


def scan_genome(tmp, records, guides, max_mm, threads=2, extra=()):
    """pack `records` [(name, seq)...], write guides, scan, return rows."""
    fa = os.path.join(tmp, "g.fa")
    cgb = os.path.join(tmp, "g.cgb")
    gtsv = os.path.join(tmp, "guides.tsv")
    out = os.path.join(tmp, "out.tsv")
    write_fasta(fa, records)
    run("pack", fa, cgb)
    write_guides(gtsv, guides)
    run("scan", cgb, gtsv, max_mm, out, threads, *extra)
    return parse_tsv(out)[1]


def rows_as_tuples(rows, keys):
    return [tuple(r[k] for k in keys) for r in rows]


SCAN_KEYS = ("guide_id", "system", "chrom", "pos", "strand", "pam",
             "mismatches", "plus_window")

# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------


def t01_pack_unpack_roundtrip():
    """pack/unpack roundtrip: lowercase, N runs, blank lines, CRLF, multi-line."""
    with tmpdir() as d:
        fa = os.path.join(d, "in.fa")
        out_fa = os.path.join(d, "out.fa")
        cgb = os.path.join(d, "in.cgb")
        a1 = "acgtACGTacgt"                 # first line of chrA
        a2 = "acgtacgtacgtNNNNNacgtNN"       # second line of chrA
        seqb = "".join("ACGTNacgtn"[i % 10] for i in range(137))
        recs = [
            ("chrA some description here", a1 + a2),
            ("chrB", seqb),
            ("chrC", "ACGT"),
            ("chrD", ""),  # empty record: allowed, must survive roundtrip
        ]
        # raw FASTA with CRLF, blank lines, and multi-line sequences
        parts = []
        parts.append(">chrA some description here")
        parts.append(a1)          # first chunk of chrA
        parts.append("")          # blank line mid-record
        parts.append(a2)
        parts.append("")          # blank line between records
        parts.append(">chrB")
        for i in range(0, len(seqb), 40):
            parts.append(seqb[i : i + 40])
        parts.append(">chrC")
        parts.append("ACGT")
        parts.append("")
        parts.append(">chrD")
        raw = "\r\n".join(parts) + "\r\n"
        with open(fa, "wb") as f:
            f.write(raw.encode())

        run("pack", fa, cgb)
        run("unpack", cgb, out_fa)
        got = read_fasta(out_fa)
        # names are the first whitespace-delimited token of the '>' line
        want = [(name.split()[0], seq.upper()) for name, seq in recs]
        assert got == want, f"roundtrip mismatch:\n got={got}\nwant={want}"


def t02_plus_strand_spcas9_exact():
    """plus-strand SPCAS9 exact site: guide+'AGG' planted -> exactly 1 site."""
    rng = random.Random(20261022)
    g = rand_seq(rng, 20)
    pos = 3000
    c1 = list(rand_seq(rng, 6000))
    c1[pos : pos + 23] = list(g + "AGG")
    records = [("c1", "".join(c1)), ("c2", rand_seq(rng, 4000))]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0)
        header, _ = parse_tsv(os.path.join(d, "out.tsv"))
        assert header == list(SCAN_KEYS), f"unexpected header {header}"
        assert len(rows) == 1, f"expected exactly 1 row, got {rows}"
        want = ("g1", "SPCAS9", "c1", str(pos + 1), "+", "AGG", "0", g)
        assert rows_as_tuples(rows, SCAN_KEYS) == [want], rows


def t03_pam_distal_mismatches():
    """2 mismatches at guide 5' end (PAM-distal): found at mm<=4, not at mm=1."""
    rng = random.Random(20261023)
    g = rand_seq(rng, 20)
    pos = 3000
    c1 = list(rand_seq(rng, 6000))
    c1[pos : pos + 23] = list(g + "AGG")
    records = [("c1", "".join(c1)), ("c2", rand_seq(rng, 4000))]
    g2 = mutate(g, [0, 1])  # guide positions 1-2 (1-based), the 5' end
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g2", "SPCAS9", g2)], 4)
        want = ("g2", "SPCAS9", "c1", str(pos + 1), "+", "AGG", "2", g)
        assert rows_as_tuples(rows, SCAN_KEYS) == [want], rows
        # not found when max_mm=1
        rows = scan_genome(d, records, [("g2", "SPCAS9", g2)], 1)
        assert rows == [], f"2-mm site must not be found at max_mm=1: {rows}"


def t04_minus_strand_worked_example():
    """minus-strand NGG site per the spec worked example: 'CC'+'A'+rc(g)."""
    rng = random.Random(20261024)
    g = rand_seq(rng, 20)
    assert g == "AAACCCGGGTTTACGTACGT" or True  # example uses any 20nt
    pos = 4000
    c1 = list(rand_seq(rng, 6000))
    c1[pos : pos + 23] = list("CCA" + revcomp(g))
    records = [("c1", "".join(c1)), ("c2", rand_seq(rng, 4000))]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0)
        want = ("g1", "SPCAS9", "c1", str(pos + 4), "-", "TGG", "0",
                revcomp(g))
        assert rows_as_tuples(rows, SCAN_KEYS) == [want], rows


def t05_nag_pam():
    """NAG PAM found by default, suppressed by --ngg-only."""
    rng = random.Random(20261025)
    g = rand_seq(rng, 20)
    pos = 2500
    c1 = list(rand_seq(rng, 5000))
    c1[pos : pos + 23] = list(g + "CAG")
    records = [("c1", "".join(c1))]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0)
        want = ("g1", "SPCAS9", "c1", str(pos + 1), "+", "CAG", "0", g)
        assert rows_as_tuples(rows, SCAN_KEYS) == [want], rows
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0,
                           extra=("--ngg-only",))
        assert rows == [], f"NAG site must be suppressed by --ngg-only: {rows}"


def t06_genome_n_in_window():
    """genome N inside window: absent at mm=0, mm=1 at max_mm=1."""
    rng = random.Random(20261026)
    g = rand_seq(rng, 20)
    pos = 3000
    c1 = list(rand_seq(rng, 6000))
    c1[pos : pos + 23] = list(g + "AGG")
    c1[pos + 10] = "N"  # 11th base of the window
    records = [("c1", "".join(c1))]
    want_window = g[:10] + "N" + g[11:]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0)
        assert rows == [], f"N-in-window must not match at mm=0: {rows}"
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 1)
        want = ("g1", "SPCAS9", "c1", str(pos + 1), "+", "AGG", "1",
                want_window)
        assert rows_as_tuples(rows, SCAN_KEYS) == [want], rows


def t07_sacas9_pams():
    """SACAS9 canonical/noncanonical PAMs on both strands; wrong PAM rejected.

    Note on the spec's test list: it says NNGRRA (last base T->A) must NOT
    be found.  Re-derived from first principles (see src/crispr_scan.c
    header): noncanonical NNGRRV = NNGRR{A,C,G} includes A, so NNGRRA IS a
    valid noncanonical PAM and must be found by default; it is only
    excluded under --canonical-only.  This test asserts both behaviors.
    A truly wrong PAM (R->T at position 4: NNGTAT... i.e. 'TTGTAT') fails
    both canonical and noncanonical and is never found.
    """
    rng = random.Random(20261027)
    s = rand_seq(rng, 21)
    rcs = revcomp(s)
    c1 = list(rand_seq(rng, 9000))
    # plus canonical NNGRRT: N N G R R T
    c1[1000 : 1000 + 27] = list(s + "TTGAAT")
    # plus noncanonical NNGRRV (last base C)
    c1[1500 : 1500 + 27] = list(s + "TTGAAC")
    # minus canonical: plus shows A T T C G A + rc(s)  -> minus PAM TCGAAT
    c1[2000 : 2000 + 27] = list("ATTCGA" + rcs)
    # minus noncanonical: plus shows G T T C G A + rc(s) -> minus PAM TCGAAC
    c1[2500 : 2500 + 27] = list("GTTCGA" + rcs)
    # wrong PAM: R->T at position 4 -> fails canonical AND noncanonical
    c1[3000 : 3000 + 27] = list(s + "TTGTAT")
    # NNGRRA (spec's example): valid NONCANONICAL only
    c1[3500 : 3500 + 27] = list(s + "TTGAAA")
    records = [("c1", "".join(c1))]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("s1", "SACAS9", s)], 0)
        want = [
            ("s1", "SACAS9", "c1", "1001", "+", "TTGAAT", "0", s),
            ("s1", "SACAS9", "c1", "1501", "+", "TTGAAC", "0", s),
            ("s1", "SACAS9", "c1", "2007", "-", "TCGAAT", "0", rcs),
            ("s1", "SACAS9", "c1", "2507", "-", "TCGAAC", "0", rcs),
            ("s1", "SACAS9", "c1", "3501", "+", "TTGAAA", "0", s),
        ]
        assert rows_as_tuples(rows, SCAN_KEYS) == want, rows
        # --canonical-only keeps NNGRRT only
        rows = scan_genome(d, records, [("s1", "SACAS9", s)], 0,
                           extra=("--canonical-only",))
        want = [
            ("s1", "SACAS9", "c1", "1001", "+", "TTGAAT", "0", s),
            ("s1", "SACAS9", "c1", "2007", "-", "TCGAAT", "0", rcs),
        ]
        assert rows_as_tuples(rows, SCAN_KEYS) == want, rows


def t08_edge_cases():
    """boundary sites, multi-chrom hits, duplicated sites."""
    rng = random.Random(20261028)
    g = rand_seq(rng, 20)
    rcs = revcomp(g)
    records = [
        ("e1", "CA" + rcs),            # minus site truncated at start: NO
        ("e2", g + "AG"),              # plus site truncated at end: NO
        ("e3", "CCA" + rcs),           # minus site fits exactly: YES (pos 4)
        ("e4", g + "AGG"),             # plus site fits exactly: YES (pos 1)
        ("e5", rand_seq(rng, 1000)),   # two copies of the site in one chrom
        ("e6", rand_seq(rng, 900)),    # one copy -> per-chrom reporting
    ]
    recs = dict(records)
    recs["e5"] = recs["e5"][:100] + g + "AGG" + recs["e5"][100:]
    recs["e5"] = recs["e5"][:500] + g + "AGG" + recs["e5"][500:]
    recs["e6"] = recs["e6"][:150] + g + "AGG" + recs["e6"][150:]
    records = [(name, recs[name]) for name, _ in records]
    with tmpdir() as d:
        rows = scan_genome(d, records, [("g1", "SPCAS9", g)], 0)
        want = [
            ("g1", "SPCAS9", "e3", "4", "-", "TGG", "0", rcs),
            ("g1", "SPCAS9", "e4", "1", "+", "AGG", "0", g),
            ("g1", "SPCAS9", "e5", "101", "+", "AGG", "0", g),
            ("g1", "SPCAS9", "e5", "501", "+", "AGG", "0", g),
            ("g1", "SPCAS9", "e6", "151", "+", "AGG", "0", g),
        ]
        assert rows_as_tuples(rows, SCAN_KEYS) == want, rows


def t09_kmap():
    """kmap: FWD + RC hits, absent k-mer, N-skip, genome N breaks match."""
    rng = random.Random(20261029)
    X = rand_seq(rng, 20)     # planted FWD in k1
    Y = rand_seq(rng, 20)     # planted as revcomp in k2
    Z = rand_seq(rng, 20)     # absent
    P = X[:10] + revcomp(X[:10])  # palindromic 20-mer
    assert P == revcomp(P)
    Xn = X[:5] + "N" + X[6:]  # genome copy of X broken by an N

    k1 = list(rand_seq(rng, 2000))
    k1[1000 : 1020] = list(X)
    k1[2000 : 2020] = list(P)
    k1[3000 : 3020] = list(Xn)
    k2 = list(rand_seq(rng, 3000))
    k2[1500 : 1520] = list(revcomp(Y))
    records = [("k1", "".join(k1)), ("k2", "".join(k2))]

    query = [
        ("qx", X),                          # FWD hit at k1:1001
        ("qy", Y),                          # RC hit at k2:1501 (kmer=rc(Y))
        ("qz", Z),                          # absent: no rows
        ("qoff", "A" * 10 + X),             # X at query_pos 11
        ("qn", X[:4] + "N" + X[5:]),        # k-mer with N: skipped
        ("qp", P),                          # palindrome: FWD and RC rows
    ]
    kmap_keys = ("query_name", "query_pos", "orient", "chrom", "chrom_pos",
                 "kmer")
    with tmpdir() as d:
        fa = os.path.join(d, "g.fa")
        cgb = os.path.join(d, "g.cgb")
        qfa = os.path.join(d, "q.fa")
        out = os.path.join(d, "k.tsv")
        write_fasta(fa, records)
        run("pack", fa, cgb)
        write_fasta(qfa, query)
        run("kmap", cgb, qfa, 20, out, 2)
        header, rows = parse_tsv(out)
        assert header == list(kmap_keys), header
        want = [
            ("qx", "1", "FWD", "k1", "1001", X),
            ("qy", "1", "RC", "k2", "1501", revcomp(Y)),
            ("qoff", "11", "FWD", "k1", "1001", X),
            ("qp", "1", "FWD", "k1", "2001", P),
            ("qp", "1", "RC", "k1", "2001", P),
        ]
        assert rows_as_tuples(rows, kmap_keys) == want, rows
        assert not any(r["query_name"] == "qz" for r in rows), "absent k-mer hit"
        assert not any(r["query_name"] == "qn" for r in rows), "N k-mer skipped"


def t10_determinism():
    """scan and kmap output byte-identical across thread counts (1 vs 6)."""
    rng = random.Random(20261030)
    g = rand_seq(rng, 20)
    rcs = revcomp(g)
    c1 = list(rand_seq(rng, 4000))
    c1[500 : 523] = list(g + "AGG")
    c1[1000 : 1023] = list("CCA" + rcs)
    c1[2000 : 2023] = list(mutate(g, [19]) + "AGG")  # 1-mm plus site
    c2 = list(rand_seq(rng, 3000))
    c2[800 : 820] = list(g)
    records = [("c1", "".join(c1)), ("c2", "".join(c2))]
    guides = [("g1", "SPCAS9", g), ("g2", "SPCAS9", mutate(g, [5, 6]))]
    query = [("q", g), ("q2", revcomp(g))]
    with tmpdir() as d:
        fa = os.path.join(d, "g.fa")
        cgb = os.path.join(d, "g.cgb")
        gtsv = os.path.join(d, "guides.tsv")
        qfa = os.path.join(d, "q.fa")
        write_fasta(fa, records)
        run("pack", fa, cgb)
        write_guides(gtsv, guides)
        write_fasta(qfa, query)
        outs = []
        for t in (1, 6):
            s_out = os.path.join(d, f"s{t}.tsv")
            k_out = os.path.join(d, f"k{t}.tsv")
            run("scan", cgb, gtsv, 2, s_out, t)
            run("kmap", cgb, qfa, 20, k_out, t)
            outs.append((s_out, k_out))
        for f1, f2 in zip(outs[0], outs[1]):
            b1 = open(f1, "rb").read()
            b2 = open(f2, "rb").read()
            assert b1 == b2, f"non-deterministic output for {f1} vs {f2}"


def t11_performance():
    """perf: 20 Mbp random genome, 200 guides, mm<=4, 6 threads < 120 s."""
    rng = random.Random(1234567)
    nch, chlen = 20, 1_000_000
    total_bp = nch * chlen
    with tmpdir() as d:
        fa = os.path.join(d, "perf.fa")
        with open(fa, "w") as f:
            for c in range(nch):
                f.write(f">perf{c}\n")
                f.write("".join(rng.choices("ACGT", k=chlen)))
                f.write("\n")
        cgb = os.path.join(d, "perf.cgb")
        run("pack", fa, cgb)
        guides = [(f"g{i:03d}", "SPCAS9", rand_seq(rng, 20))
                  for i in range(200)]
        gtsv = os.path.join(d, "guides.tsv")
        out = os.path.join(d, "perf.tsv")
        write_guides(gtsv, guides)
        t0 = time.monotonic()
        run("scan", cgb, gtsv, 4, out, 6)
        dt = time.monotonic() - t0
        rows = parse_tsv(out)[1]
        extrap = dt * GRCH38_BP / total_bp
        print(f"\n  PERF: scan {total_bp/1e6:.0f} Mbp x 200 guides "
              f"(max_mm=4, both strands+PAMs) 6 threads: {dt:.2f} s; "
              f"hits={len(rows)}")
        print(f"  PERF: extrapolated GRCh38 3.1 Gbp: {extrap/60:.1f} min "
              f"(spec target < 30 min)")
        assert dt < 120, f"scan took {dt:.1f}s (budget 120s)"
        assert len(rows) > 0, "expected at least one off-target hit"


def t12_stat():
    """stat output: n_chrom, total_bases and per-chrom table."""
    rng = random.Random(20261031)
    lens = [1234, 567]
    records = [("chrA", rand_seq(rng, lens[0])), ("chrB", rand_seq(rng, lens[1]))]
    with tmpdir() as d:
        fa = os.path.join(d, "g.fa")
        cgb = os.path.join(d, "g.cgb")
        write_fasta(fa, records)
        run("pack", fa, cgb)
        r = run("stat", cgb)
        lines = r.stdout.strip().splitlines()
        assert lines[0] == f"n_chrom\t2", lines
        assert lines[1] == f"total_bases\t{sum(lens)}", lines
        assert lines[2] == f"chrA\t{lens[0]}", lines
        assert lines[3] == f"chrB\t{lens[1]}", lines


def t13_cli_errors():
    """error handling: bad inputs exit non-zero with a stderr message."""
    with tmpdir() as d:
        fa = os.path.join(d, "ok.fa")
        cgb = os.path.join(d, "ok.cgb")
        write_fasta(fa, [("c", "ACGT" * 100)])
        run("pack", fa, cgb)

        # --help exits 0 and documents all modes
        r = run("--help")
        assert r.returncode == 0 and "pack" in r.stdout and "scan" in r.stdout
        assert "kmap" in r.stdout and "stat" in r.stdout and "unpack" in r.stdout

        # bad guide characters (exactly 20 nt so the char check is what fires)
        g1 = os.path.join(d, "gN.tsv")
        write_guides(g1, [("bad", "SPCAS9", "ACGTNACGTACGTACGTACG")])
        r = run("scan", cgb, g1, 3, os.path.join(d, "o.tsv"), expect_ok=False)
        assert r.returncode != 0 and "invalid character" in r.stderr, r.stderr

        # wrong SPCAS9 guide length
        g2 = os.path.join(d, "gL.tsv")
        write_guides(g2, [("bad", "SPCAS9", "ACGT" * 5 + "A")])
        r = run("scan", cgb, g2, 3, os.path.join(d, "o.tsv"), expect_ok=False)
        assert r.returncode != 0 and "20" in r.stderr, r.stderr

        # unknown system
        g3 = os.path.join(d, "gS.tsv")
        write_guides(g3, [("bad", "CPF1", "ACGT" * 5)])
        r = run("scan", cgb, g3, 3, os.path.join(d, "o.tsv"), expect_ok=False)
        assert r.returncode != 0 and "SPCAS9" in r.stderr, r.stderr

        # bad magic
        bad = os.path.join(d, "bad.cgb")
        with open(bad, "wb") as f:
            f.write(b"this is not a cgb file at all")
        r = run("stat", bad, expect_ok=False)
        assert r.returncode != 0 and "magic" in r.stderr, r.stderr

        # truncated cgb (chop the valid one)
        trunc = os.path.join(d, "trunc.cgb")
        blob = open(cgb, "rb").read()
        with open(trunc, "wb") as f:
            f.write(blob[: len(blob) // 2])
        r = run("stat", trunc, expect_ok=False)
        assert r.returncode != 0 and "truncat" in r.stderr, r.stderr

        # k out of range
        r = run("kmap", cgb, fa, 32, os.path.join(d, "o.tsv"),
                expect_ok=False)
        assert r.returncode != 0, r.stderr
        r = run("kmap", cgb, fa, 0, os.path.join(d, "o.tsv"),
                expect_ok=False)
        assert r.returncode != 0, r.stderr

        # missing / bad arguments
        r = run("scan", cgb, g1, expect_ok=False)
        assert r.returncode != 0, r.stderr
        r = run("scan", cgb, g1, "notanumber", os.path.join(d, "o.tsv"),
                expect_ok=False)
        assert r.returncode != 0, r.stderr
        r = run("nosuchmode", expect_ok=False)
        assert r.returncode != 0, r.stderr


# ----------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------

TESTS = [t01_pack_unpack_roundtrip, t02_plus_strand_spcas9_exact,
        t03_pam_distal_mismatches, t04_minus_strand_worked_example,
        t05_nag_pam, t06_genome_n_in_window, t07_sacas9_pams,
        t08_edge_cases, t09_kmap, t10_determinism, t11_performance,
        t12_stat, t13_cli_errors]


def main():
    ensure_bin()
    failed = 0
    for t in TESTS:
        try:
            t()
            print(f"PASS {t.__name__}: "
                  f"{(t.__doc__ or '').strip().splitlines()[0]}")
        except Exception as e:
            failed += 1
            print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
    print(f"\n{len(TESTS) - failed}/{len(TESTS)} tests passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()