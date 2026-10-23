#!/usr/bin/env python3
"""Render the paper's Sources block mechanically from the citation ledger.
Extracts \\src{n} ids from paper/*.tex, maps to ledger URLs, emits
paper/sources_block.tex in evidence style (quotes where captured).
Run: /usr/bin/python3 src/make_sources_block.py
"""
import json
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(HERE, "paper")
LEDGER = os.path.join(HERE, "docs", "citations", "ledger.json")


GREEK = {"\u03b1": r"$\alpha$", "\u03b2": r"$\beta$", "\u03b3": r"$\gamma$",
        "\u03b4": r"$\delta$", "\u03b5": r"$\epsilon$", "\u03b6": r"$\zeta$",
        "\u03b7": r"$\eta$", "\u03b8": r"$\theta$", "\u03b9": r"i",
        "\u03ba": r"$\kappa$", "\u03bb": r"$\lambda$", "\u03bc": r"$\mu$",
        "\u03c3": r"$\sigma$", "\u03c4": r"$\tau$", "\u03c7": r"$\chi$",
        "\u03a9": r"$\Omega$", "\u00b1": r"$\pm$", "\u00d7": r"$\times$",
        "\u2248": r"$\approx$", "\u2264": r"$\le$", "\u2265": r"$\ge$",
        "\u00b0": r"$^\circ$", "\u00e9": r"\'e", "\u00e8": r"\'e",
        "\u00fc": r"\"u", "\u00f6": r"\"o", "\u00e4": r"\"a"}


def latex_escape(s):
    repl = [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
            ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
            ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}"),
            ("\u2018", "`"), ("\u2019", "'"), ("\u201c", "``"), ("\u201d", "''"),
            ("\u2014", "---"), ("\u2013", "--"), ("\u2026", "..."),
            ("\u2032", "$'$"), ("\u00a0", "~")]
    out = []
    for ch in s:
        if ch in GREEK:
            out.append(GREEK[ch])
            continue
        for a, b in repl:
            if ch == a:
                out.append(b)
                break
        else:
            if ord(ch) > 126:
                out.append(f"[{ch}]")  # flag unexpected unicode visibly
            else:
                out.append(ch)
    return "".join(out)


def main():
    tex = ""
    for fn in os.listdir(PAPER):
        if fn.endswith(".tex") and fn != "sources_block.tex":
            tex += open(os.path.join(PAPER, fn)).read()
    cited = sorted({int(n) for n in re.findall(r"\\src\{(\d+)\}", tex)})
    ledger = json.load(open(LEDGER))
    srcs = ledger.get("sources", ledger)
    by_id = {s["id"]: s for s in srcs}
    with open(os.path.join(PAPER, "sources_block.tex"), "w") as fh:
        fh.write("\\section*{Sources}\n\\small\n")
        fh.write("Entries below were registered at retrieval time and verified "
                 "against the fetched page text; quoted passages are verbatim "
                 "substrings of the saved page files in "
                 "\\texttt{docs/citations/pages/}. Sources without a quoted "
                 "passage are data repositories or pages that yielded no clean "
                 "quotable sentence; they are cited from their fetched "
                 "content. Very long machine-extracted passages (over 200 "
                 "characters, e.g.\\ raw sequence translations) are omitted "
                 "from the rendered quotes for typography; they remain in the "
                 "ledger file.\n\n\\begin{enumerate}\n")
        n_q = 0
        for i in cited:
            s = by_id.get(i)
            if s is None:
                raise SystemExit(f"src{{{i}}} not in ledger!")
            title = latex_escape(s.get("title", "") or "")
            url = s.get("url", "")
            quotes = s.get("quotes", [])
            fh.write(f"\\item [{i}] {title} \\\\ \\url{{{url}}}\n")
            for q in quotes:
                qt = (q.get("text") or "").strip()
                if len(qt) > 200 or max(len(w) for w in qt.split()) > 38:
                    continue  # machine-extracted noise; kept in ledger only
                qt = latex_escape(qt)
                if qt:
                    fh.write(f"  \\\\ {{\\sloppy \\emph{{Quote:}} ``{qt}''}}\n")
                    n_q += 1
            fh.write("\n")
        fh.write("\\end{enumerate}\n")
    print(f"sources_block.tex written: {len(cited)} cited sources, "
          f"{n_q} verbatim quotes")


if __name__ == "__main__":
    main()