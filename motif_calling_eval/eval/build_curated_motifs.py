#!/usr/bin/env python3
"""Build curated REBASE motif TSV: motifs in >=3 genomes, known methylation
type, median percent_detection >= 50, dropping the multi_n category.

Output columns: motif, methylated_base, methylation_type, n_genomes, median_pct, category
"""
from __future__ import annotations
import csv, statistics, argparse, sys
from collections import defaultdict, Counter
from pathlib import Path

DEFAULT_REBASE = "/home/alex/Documents/Restriction_Methylation/rebase/REBASE_Apr21_parsed.csv"

def split_n_runs(motif: str):
    parts, cur, in_n = [], [], False
    for c in motif:
        if c == "N":
            if not in_n:
                parts.append(["".join(cur), 0]); cur = []; in_n = True
            parts[-1][1] += 1
        else:
            in_n = False; cur.append(c)
    parts.append(["".join(cur), 0])
    return parts

def categorize(motif: str) -> str:
    if "N" not in motif:
        return "solid_degen" if any(c not in "ACGT" for c in motif) else "solid_simple"
    parts = split_n_runs(motif)
    segs = [(s, n) for s, n in parts if s != "" or n > 0]
    n_runs = [n for _, n in segs[:-1] if n > 0]
    halves = [s for s, _ in segs if s != ""]
    if len(n_runs) >= 2 or len(halves) != 2:
        return "multi_n"
    sp = n_runs[0]
    a, b = len(halves[0]), len(halves[1])
    if sp >= 9: return "bip_sp9plus"
    if 1 <= sp <= 3:
        return "bip_sp1_3_seedable" if max(a,b) >= 4 else "bip_sp1_3_noseed"
    return "bip_sp4_8_good" if min(a,b) >= 3 else "bip_sp4_8_shorthalf"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebase", default=DEFAULT_REBASE)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-genomes", type=int, default=3)
    ap.add_argument("--min-median-pct", type=float, default=50.0)
    ap.add_argument("--drop-categories", nargs="*", default=["multi_n"])
    args = ap.parse_args()

    motif_rows = defaultdict(list)
    with open(args.rebase) as f:
        for row in csv.DictReader(f):
            m = (row.get("motif") or "").strip().upper()
            gb = (row.get("genbank") or "").strip()
            mt = (row.get("methylation_type") or "").strip()
            mb = row.get("methylated_base", "")
            try:
                mb_int = int(float(mb))
                pct = float(row.get("percent_detection") or "")
            except (ValueError, TypeError):
                continue
            if not m or not gb: continue
            if not all(c in "ACGTRYWSKMBDHVN" for c in m): continue
            if mt in ("", "Unknown"): continue
            motif_rows[m].append((gb, mt, mb_int, pct))

    drop = set(args.drop_categories)
    out_rows = []
    for m, rs in motif_rows.items():
        gbs = {r[0] for r in rs}
        if len(gbs) < args.min_genomes: continue
        med = statistics.median(r[3] for r in rs)
        if med < args.min_median_pct: continue
        cat = categorize(m)
        if cat in drop: continue
        # most common methylation type and methylated_base for this motif
        mt_mode = Counter(r[1] for r in rs).most_common(1)[0][0]
        # mode of mb among rows that match the modal methylation type
        mb_mode = Counter(r[2] for r in rs if r[1] == mt_mode).most_common(1)[0][0]
        out_rows.append((m, mb_mode, mt_mode, len(gbs), med, cat))

    out_rows.sort(key=lambda x: (-x[3], x[0]))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write("motif\tmethylated_base\tmethylation_type\tn_genomes\tmedian_pct\tcategory\n")
        for r in out_rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.1f}\t{r[5]}\n")
    cat_counts = Counter(r[5] for r in out_rows)
    print(f"Wrote {len(out_rows)} motifs to {args.out}", file=sys.stderr)
    for k, v in sorted(cat_counts.items()):
        print(f"  {k:28s} {v}", file=sys.stderr)

if __name__ == "__main__":
    main()
