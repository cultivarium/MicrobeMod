#!/usr/bin/env python3
"""Production-faithful evaluation harness for microbe_motif.

Differences from the old eval_motif_caller.py (which were wrong):
  1. PRODUCTION POST-PROCESSING. MicrobeMod does not use find_motifs' raw IUPAC
     output. It writes streme.xml (per-position PWM regenerated from the IUPAC
     char via _iupac_pwm_row) and parses it with assign_motifs(), which:
        - drops motifs with e-value >= MIN_EVALUE (0.1)
        - re-masks 3-letter codes (B/D/H/V) -> N   (1/2-letter codes survive,
          since their regenerated PWM mass is 0.97 > MOTIF_FREQ_CUTOFF=0.8)
        - strips leading/trailing N
     We replicate exactly that transform on each called motif.
  2. ALL MOD TYPES PER GENOME. Production (microbemod.py:586) loops over every
     methylation type, building pos.fasta from only that mod's sites and calling
     find_motifs once per type. We do the same and UNION the called motifs, then
     score one P/R/F1 over ALL truth motifs of the genome.
  3. SAVE MOTIF OUTPUTS. Per replicate we dump truth + raw + post-processed
     called motifs (per mod type) to motifs/<rep>.json.
  4. DETERMINISTIC. sorted() mod ordering + PYTHONHASHSEED (set by the launcher).

Which commit is exercised is chosen by MM_ROOT (sys.path).
"""
from __future__ import annotations
import argparse, json, os, sys, time
from collections import defaultdict
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

# Which MicrobeMod to import (commit under test). Set via MM_ROOT env so it
# propagates to forkserver/spawn workers; falls back to this repo.
MM_ROOT = os.environ.get("MM_ROOT", str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, MM_ROOT)
from MicrobeMod import microbe_motif as mm  # noqa: E402

MOD_CODE_INV = {"6mA": "a", "5mC": "m", "4mC": "21839", "5hmC": "h"}
MATCH_MAX_H = 2
# Production constants (microbemod.py)
MIN_EVALUE = 0.1
THREE_LETTER = set("BDHV")


def production_transform(iupac: str, evalue: float) -> str | None:
    """Replicate microbemod.assign_motifs() motif-string handling. Returns the
    production motif string, or None if it would be dropped."""
    if evalue >= MIN_EVALUE:
        return None
    masked = "".join("N" if c in THREE_LETTER else c for c in iupac)
    stripped = masked.strip("N")
    return stripped or None


def f1_of(called: list[str], truth: list[str]) -> dict:
    cm = [False] * len(called)
    tm = [False] * len(truth)
    for ci, c in enumerate(called):
        for ti, t in enumerate(truth):
            if mm.motifs_alignable(c, t, MATCH_MAX_H):
                cm[ci] = True
                tm[ti] = True
    tp = sum(cm)
    fp = len(called) - tp
    fn = len(truth) - sum(tm)
    p = tp / (tp + fp) if (tp + fp) else 0.0
    r = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"tp": tp, "fp": fp, "fn": fn,
            "precision": round(p, 4), "recall": round(r, 4), "f1": round(f1, 4)}


def evaluate_one(combo_dir, genome_path, truth, bg, motif_dir):
    bed_name = truth.get("bed")
    bed_path = os.path.join(combo_dir, bed_name)
    if not os.path.exists(bed_path):
        return None
    truth_motifs = [m["motif"] for m in truth.get("motifs", [])]
    k = len(truth_motifs)
    # Mod types present, deterministic order.
    mod_types = sorted({m.get("methylation_type", "6mA") for m in truth.get("motifs", [])})

    per_mod = []
    called_all = []
    t_find = 0.0  # wall time spent inside find_motifs (the work under test)
    for mt in mod_types:
        mod_code = MOD_CODE_INV.get(mt, mt)
        try:
            pos, neg = mm._build_pos_neg_from_bed(
                bed_path, genome_path, mod_code,
                percent_cutoff=0.0, min_coverage=1, n_control=10000, seed=13)
        except Exception as e:
            per_mod.append({"mod": mt, "error": f"build:{e}"})
            continue
        if pos.shape[0] < mm.MIN_MOTIF_SITES:
            per_mod.append({"mod": mt, "n_pos": int(pos.shape[0]),
                            "raw": [], "prod": []})
            continue
        try:
            _t0 = time.time()
            motifs = mm.find_motifs(pos, neg, bg=bg)
            t_find += time.time() - _t0
        except Exception as e:
            per_mod.append({"mod": mt, "n_pos": int(pos.shape[0]),
                            "error": f"find:{e}"})
            continue
        raw = [m.iupac for m in motifs]
        prod = []
        for m in motifs:
            t = production_transform(m.iupac, float(m.evalue))
            if t:
                prod.append(t)
        called_all.extend(prod)
        per_mod.append({"mod": mt, "n_pos": int(pos.shape[0]),
                        "raw": raw, "prod": prod})

    metrics = f1_of(called_all, truth_motifs)
    rec = {"bed": bed_name, "k": k, "mod_types": mod_types,
           "n_called": len(called_all), **metrics,
           "time_s": round(t_find, 3), "n_mod_runs": len(mod_types),
           "called": called_all, "truth": truth_motifs, "per_mod": per_mod}
    # Save per-rep motif outputs.
    (Path(motif_dir) / (bed_name.replace(".bed", ".json"))).write_text(
        json.dumps(rec, indent=2))
    return rec


def k_group(k):
    return "1-2" if k <= 2 else "3-5" if k <= 5 else "6-10" if k <= 10 else ">10"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combo-dir", default=str(Path(__file__).resolve().parents[1] / "eval_output/combo"))
    ap.add_argument("--genome", default=str(Path(__file__).resolve().parents[2] / "tests/test_data/EcoliCVM05_GCF_000005845.2_ASM584v2_genomic.fna"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--label", default="run")
    ap.add_argument("--max-workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    combo_dir = Path(args.combo_dir)
    genome_path = str(Path(args.genome).resolve())
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    motif_dir = out_dir / "motifs"; motif_dir.mkdir(exist_ok=True)

    truth_files = sorted(combo_dir.glob("rep*.truth.json"))
    if args.limit:
        truth_files = truth_files[:args.limit]
    print(f"[{args.label}] MM_ROOT={MM_ROOT}", file=sys.stderr)
    print(f"[{args.label}] {len(truth_files)} reps, {args.max_workers} workers",
          file=sys.stderr)

    bg = mm.compute_genome_bg(genome_path)
    results, errors = [], []
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=args.max_workers) as exc:
        futs = {}
        for tf in truth_files:
            truth = json.loads(tf.read_text())
            if "bed" not in truth:
                truth["bed"] = tf.name.replace(".truth.json", ".bed")
            futs[exc.submit(evaluate_one, str(combo_dir), genome_path, truth,
                            bg, str(motif_dir))] = tf.name
        for fut in as_completed(futs):
            try:
                r = fut.result()
                if r:
                    results.append(r)
            except Exception as e:
                errors.append({"file": futs[fut], "error": str(e)})
    total_wall_s = time.time() - t_start

    groups = defaultdict(list)
    for r in results:
        groups[k_group(r["k"])].append(r)
    print(f"\n=== [{args.label}] {len(results)} reps, {len(errors)} errors ===")
    print(f"{'group':<8}{'n':>4}{'meanF1':>9}{'medF1':>9}{'meanP':>9}{'meanR':>9}")
    summary = {}
    for g in ["1-2", "3-5", "6-10", ">10", "ALL"]:
        rs = results if g == "ALL" else groups.get(g, [])
        if not rs:
            continue
        f1s = [r["f1"] for r in rs]
        summary[g] = {"n": len(rs),
                      "mean_f1": round(float(np.mean(f1s)), 4),
                      "median_f1": round(float(np.median(f1s)), 4),
                      "mean_precision": round(float(np.mean([r["precision"] for r in rs])), 4),
                      "mean_recall": round(float(np.mean([r["recall"] for r in rs])), 4)}
        s = summary[g]
        print(f"{g:<8}{s['n']:>4}{s['mean_f1']:>9.4f}{s['median_f1']:>9.4f}"
              f"{s['mean_precision']:>9.4f}{s['mean_recall']:>9.4f}")

    sum_find = sum(r.get("time_s", 0.0) for r in results)
    sum_runs = sum(r.get("n_mod_runs", 0) for r in results)
    summary["_meta"] = {
        "label": args.label, "mm_root": MM_ROOT,
        "n_reps": len(results), "max_workers": args.max_workers,
        "total_wall_s": round(total_wall_s, 1),
        "sum_find_motifs_s": round(sum_find, 1),
        "n_find_motifs_calls": sum_runs,
        "mean_find_s_per_call": round(sum_find / sum_runs, 3) if sum_runs else 0.0,
        "mean_wall_s_per_rep": round(total_wall_s / len(results), 2) if results else 0.0,
    }
    print(f"[{args.label}] wall={total_wall_s:.1f}s  "
          f"sum_find_motifs={sum_find:.1f}s over {sum_runs} calls  "
          f"({summary['_meta']['mean_find_s_per_call']}s/call)", file=sys.stderr)

    import csv
    with open(out_dir / "eval_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["bed", "k", "group", "n_called",
                           "tp", "fp", "fn", "precision", "recall", "f1",
                           "time_s", "n_mod_runs"],
                           extrasaction="ignore")
        w.writeheader()
        for r in results:
            r["group"] = k_group(r["k"])
            w.writerow(r)
    (out_dir / "eval_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[{args.label}] wrote {out_dir}/ (motifs/, eval_results.csv, eval_summary.json)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
