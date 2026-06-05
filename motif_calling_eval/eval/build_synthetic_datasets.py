#!/usr/bin/env python3
"""Synthetic methylation dataset generator for motif-caller evaluation.

For a list of motifs (with methylated_base and methylation_type), scans both
strands of a reference genome and produces modkit-compatible bedMethyl files
suitable for the MicrobeMod motif-finding pipeline.

Three dataset modes:
  single    one motif per BED file, every motif in the list
  combo     k motifs per BED file (k uniform in [k_min, k_max]); n replicates
  partial   like combo, but each motif methylated at p in [p_min, p_max]

Output layout (under <outdir>):
  single/
    <motif>__mpos<mb>.bed
    <motif>__mpos<mb>.truth.json
    ...
  combo/
    rep00001.bed
    rep00001.truth.json
    ...
  partial/
    rep00001.bed
    rep00001.truth.json
    ...

Each .truth.json describes the motifs included, methylation rates, methylated
positions count, and the contributing motif(s) for each BED row.
"""
from __future__ import annotations
import argparse, csv, json, random, re, sys
from collections import defaultdict
from pathlib import Path

# -- IUPAC --------------------------------------------------------------------
IUPAC = {
    "A":"A","C":"C","G":"G","T":"T",
    "R":"AG","Y":"CT","W":"AT","S":"CG","K":"GT","M":"AC",
    "B":"CGT","D":"AGT","H":"ACT","V":"ACG","N":"ACGT",
}
RC_TBL = str.maketrans("ACGTRYWSKMBDHVNacgtrywskmbdhvn",
                       "TGCAYRWSMKVHDBNtgcayrwsmkvhdbn")
def revcomp(s: str) -> str: return s.translate(RC_TBL)[::-1]
def iupac_re(motif: str) -> re.Pattern:
    return re.compile("".join("[" + IUPAC[c] + "]" for c in motif))

# modkit modification codes (subset; extend as needed)
MOD_CODE = {"6mA": "a", "5mC": "m", "4mC": "21839", "5hmC": "h"}

# -- Genome -------------------------------------------------------------------
def load_genome(path: Path):
    contigs = {}
    name, buf = None, []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            contigs[name] = "".join(buf).upper()
    return contigs

# -- Motif scanning -----------------------------------------------------------
def find_motif_sites(motif: str, mpos1: int, contigs: dict[str, str]):
    """Return list of (contig, methylated_pos_0idx, strand) where motif occurs.

    + strand: motif appears in contig sequence at position p; methylated base is
              at p + (mpos1 - 1).
    - strand: revcomp(motif) appears in contig at position p; methylated base
              on - strand corresponds to + strand position
              p + (len(motif) - 1) - (mpos1 - 1) = p + len(motif) - mpos1.
    """
    fwd_re = iupac_re(motif)
    rev_re = iupac_re(revcomp(motif))
    sites = []
    L = len(motif)
    for cname, seq in contigs.items():
        for m in fwd_re.finditer(seq):
            sites.append((cname, m.start() + (mpos1 - 1), "+", "fwd"))
        for m in rev_re.finditer(seq):
            sites.append((cname, m.start() + (L - mpos1), "-", "rev"))
    return sites

# -- BED writing --------------------------------------------------------------
BED_TEMPLATE = (
    "{chrom}\t{start}\t{end}\t{mod}\t{score}\t{strand}\t"
    "{start}\t{end}\t0,0,0\t{cov}\t{frac:.2f}\t{nmod}\t{ncan}\t0\t0\t{nfail}\t0\t0\n"
)

def write_bed(rows, path: Path, coverage=30, default_pct=1.0):
    """rows: iterable of dicts with keys chrom, pos, strand, mod, pct (0..1)."""
    with path.open("w") as f:
        for r in rows:
            pct = r.get("pct", default_pct)
            nmod = round(coverage * pct)
            ncan = coverage - nmod
            f.write(BED_TEMPLATE.format(
                chrom=r["chrom"], start=r["pos"], end=r["pos"]+1,
                mod=r["mod"], score=coverage, strand=r["strand"],
                cov=coverage, frac=100.0 * pct,
                nmod=nmod, ncan=ncan, nfail=0,
            ))

# -- Dataset modes ------------------------------------------------------------
def load_motifs(tsv: Path):
    out = []
    with tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out.append({
                "motif": row["motif"],
                "mpos": int(row["methylated_base"]),
                "mtype": row["methylation_type"],
                "category": row.get("category", ""),
            })
    return out

def collect_sites_for(motifs, contigs, cache):
    """Resolve all motif sites once; return cache[motif_key] = list of sites."""
    for m in motifs:
        key = (m["motif"], m["mpos"])
        if key not in cache:
            cache[key] = find_motif_sites(m["motif"], m["mpos"], contigs)
    return cache

def site_dict(site, mod_code, pct):
    cname, pos, strand, _ = site
    return {"chrom": cname, "pos": pos, "strand": strand, "mod": mod_code, "pct": pct}

def merge_collisions(rows):
    """If two motifs methylate the same (chrom,pos,strand,mod), keep the higher pct."""
    by_key = {}
    for r in rows:
        k = (r["chrom"], r["pos"], r["strand"], r["mod"])
        if k not in by_key or r["pct"] > by_key[k]["pct"]:
            by_key[k] = r
    return list(by_key.values())

def slug(motif): return motif.replace("/", "_")

def gen_single(motifs, contigs, outdir: Path, coverage: int, cache):
    outdir.mkdir(parents=True, exist_ok=True)
    truth_index = []
    collect_sites_for(motifs, contigs, cache)
    for m in motifs:
        key = (m["motif"], m["mpos"])
        sites = cache[key]
        mod_code = MOD_CODE.get(m["mtype"], m["mtype"])
        rows = [site_dict(s, mod_code, 1.0) for s in sites]
        if not rows:
            continue
        stem = f"{slug(m['motif'])}__mpos{m['mpos']}"
        bed = outdir / f"{stem}.bed"
        write_bed(rows, bed, coverage=coverage, default_pct=1.0)
        truth = {
            "mode": "single",
            "bed": bed.name,
            "motifs": [{
                "motif": m["motif"], "methylated_base": m["mpos"],
                "methylation_type": m["mtype"], "mod_code": mod_code,
                "category": m["category"],
                "n_sites": len(sites), "methylation_rate": 1.0,
            }],
            "total_rows": len(rows),
        }
        (outdir / f"{stem}.truth.json").write_text(json.dumps(truth, indent=2))
        truth_index.append({"bed": bed.name, "motif": m["motif"], "mpos": m["mpos"],
                            "n_sites": len(sites)})
    (outdir / "_index.json").write_text(json.dumps(truth_index, indent=2))
    print(f"[single] wrote {len(truth_index)} BEDs to {outdir}", file=sys.stderr)

def gen_combo(motifs, contigs, outdir: Path, n: int, k_range, coverage: int,
              cache, partial_pct_range=None, seed=0, prefix="rep"):
    """Generate n replicates; each picks k uniform in k_range motifs.
    If partial_pct_range is None, each motif methylated 100%.
    Otherwise, each motif drawn p ~ U(partial_pct_range), and a Bernoulli(p)
    sample of its sites is included (the rest are simply omitted from the BED)."""
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    collect_sites_for(motifs, contigs, cache)
    truth_index = []
    width = len(str(n))
    for i in range(1, n+1):
        k = rng.randint(k_range[0], k_range[1])
        chosen = rng.sample(motifs, k)
        rows = []
        truth_motifs = []
        for m in chosen:
            key = (m["motif"], m["mpos"])
            sites = cache[key]
            if partial_pct_range is None:
                pct = 1.0
                kept = sites
            else:
                pct = rng.uniform(partial_pct_range[0], partial_pct_range[1])
                # Bernoulli sample which sites are methylated; record per-site pct
                kept = [s for s in sites if rng.random() < pct]
            mod_code = MOD_CODE.get(m["mtype"], m["mtype"])
            rows.extend(site_dict(s, mod_code, pct) for s in kept)
            truth_motifs.append({
                "motif": m["motif"], "methylated_base": m["mpos"],
                "methylation_type": m["mtype"], "mod_code": mod_code,
                "category": m["category"],
                "n_sites_total": len(sites),
                "n_sites_methylated": len(kept),
                "methylation_rate": round(pct, 4),
            })
        rows = merge_collisions(rows)
        rep_id = f"{prefix}{i:0{width}d}"
        bed = outdir / f"{rep_id}.bed"
        write_bed(rows, bed, coverage=coverage)
        truth = {
            "mode": "combo" if partial_pct_range is None else "partial",
            "rep_id": rep_id, "bed": bed.name,
            "k": k, "seed": seed, "rng_index": i,
            "partial_pct_range": partial_pct_range,
            "motifs": truth_motifs,
            "total_rows": len(rows),
        }
        (outdir / f"{rep_id}.truth.json").write_text(json.dumps(truth, indent=2))
        truth_index.append({"bed": bed.name, "k": k,
                            "motifs": [m["motif"] for m in chosen],
                            "n_rows": len(rows)})
    (outdir / "_index.json").write_text(json.dumps(truth_index, indent=2))
    print(f"[{outdir.name}] wrote {n} replicates (k in {k_range}, partial={partial_pct_range})",
          file=sys.stderr)

# -- CLI ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--motifs", required=True, help="TSV with columns motif, methylated_base, methylation_type, category")
    ap.add_argument("--genome", required=True, help="Reference FASTA (single or multi-contig)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mode", choices=["single","combo","partial","all"], default="all")
    ap.add_argument("--coverage", type=int, default=30)
    ap.add_argument("--n-combo", type=int, default=1000)
    ap.add_argument("--n-partial", type=int, default=500)
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=20)
    ap.add_argument("--p-min", type=float, default=0.5)
    ap.add_argument("--p-max", type=float, default=0.9)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    motifs = load_motifs(Path(args.motifs))
    print(f"Loaded {len(motifs)} motifs", file=sys.stderr)
    contigs = load_genome(Path(args.genome))
    total_bp = sum(len(s) for s in contigs.values())
    print(f"Loaded {len(contigs)} contig(s), {total_bp:,} bp total", file=sys.stderr)

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    cache: dict[tuple[str,int], list] = {}

    if args.mode in ("single","all"):
        gen_single(motifs, contigs, out / "single", args.coverage, cache)
    if args.mode in ("combo","all"):
        gen_combo(motifs, contigs, out / "combo",
                  n=args.n_combo, k_range=(args.k_min, args.k_max),
                  coverage=args.coverage, cache=cache,
                  partial_pct_range=None, seed=args.seed, prefix="rep")
    if args.mode in ("partial","all"):
        gen_combo(motifs, contigs, out / "partial",
                  n=args.n_partial, k_range=(args.k_min, args.k_max),
                  coverage=args.coverage, cache=cache,
                  partial_pct_range=(args.p_min, args.p_max),
                  seed=args.seed + 1, prefix="rep")

if __name__ == "__main__":
    main()
