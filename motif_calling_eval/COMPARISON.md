# Commit comparison: motif-similarity F1 of `microbe_motif` (production-faithful)

**Question.** How does commit `4cd2529` ("reproducibility, robustness, tests"),
which replaced the set-level optimizer's wall-clock budget with a fixed pass cap,
change `microbe_motif`'s motif-similarity F1 — as MicrobeMod actually runs it?

## Harness (corrected)

`motif_calling_eval/eval/eval_harness.py` reproduces the production path:

1. **Production post-processing.** MicrobeMod does not score `find_motifs`' raw
   IUPAC output. It writes streme.xml (per-position PWM regenerated from the
   IUPAC char) and parses it with `assign_motifs()`, which drops e-value ≥ 0.1,
   re-masks 3-letter codes B/D/H/V → N (1/2-letter codes survive), and **strips
   leading/trailing N**. The harness applies exactly this transform before
   matching. (This is why the earlier "N-padding" finding was an artifact —
   production strips those N's.)
2. **All mod-types per genome.** Production (`microbemod.py:586`) loops over each
   methylation type, building pos.fasta from only that mod's sites and calling
   `find_motifs` once per type. The harness does the same and **unions** the
   called motifs, scoring one P/R/F1 over all truth motifs. (21/40 reps here are
   multi-mod. The old eval picked ONE mod-type via `list(mod_codes)[0]` — hash
   nondeterministic — and scored against all truths, capping recall. Example:
   rep002 went 0.33 → 1.00 once fixed.)
3. **Motif outputs + timing saved** per replicate (`<run>/motifs/<rep>.json`:
   truth, raw, post-processed `prod`; plus per-rep `time_s` and per-run
   `total_wall_s` in `eval_summary.json`).
4. **Deterministic** (sorted mod order, PYTHONHASHSEED=0, seeded sampling).

Three optimizer strategies, same 40 replicates, paired:

| run | strategy | commit |
|---|---|---|
| A | no optimization (`opt_max_passes=0`) | HEAD |
| B | old commit (`opt_max_seconds=10`, wall-clock) | HEAD~1 (7cd5b45) |
| C | new commit (`opt_max_passes=5`) | HEAD (4cd2529) |

## Result — F1 is essentially identical; runtime is not

| config | ALL F1 | median | meanP | meanR | wall | s / find_motifs call |
|---|---|---|---|---|---|---|
| A no-opt        | 0.9276 | 1.000 | 0.9422 | 0.9158 | 248 s | 12.7 |
| B old (10 s)    | 0.9317 | 1.000 | 0.9484 | 0.9183 | 312 s | 19.1 |
| D new (2 pass)  | 0.9300 | 1.000 | 0.9514 | 0.9136 | 925 s | 87.1 |
| C new (5 pass)  | 0.9302 | 1.000 | 0.9489 | 0.9153 | **2215 s** | **207.5** |

`passes=2` vs `passes=5`: F1 −0.0002 (35/40 ties), **2.4× faster** (925 s vs
2215 s) — confirming the 5-pass default is pure wasted runtime.

F1 by motif-count group:

| group | no-opt | old (10s) | new (5p) |
|------|--------|-----------|----------|
| 1–2  | 0.9167 | 0.9167 | 0.9167 |
| 3–5  | 0.9437 | 0.9437 | 0.9437 |
| 6–10 | 0.9427 | 0.9624 | 0.9536 |
| >10  | 0.9102 | 0.9048 | 0.9081 |

Paired (n=40): old−noopt **+0.0040**, new−noopt **+0.0026**, **new−old −0.0015**
(new better in 3 reps, worse in 5, tied in 32).

## Conclusion

- **The commit does not meaningfully change F1.** New vs old = −0.0015 (32/40
  ties). All three strategies sit within 0.004 of each other.
- **The set-level optimizer barely matters at all.** No-optimization F1 (0.9276)
  is within 0.004 of both optimized runs. Its only visible benefit is a small
  bump in the 6–10-motif group (+0.02 over no-opt); elsewhere it's flat or
  slightly negative.
- **The commit's real effect is runtime, in the wrong direction.** The 5-pass
  cap runs the optimizer to near-convergence on every dense replicate, making it
  **~7× slower wall-clock (2215 s vs 312 s)** and **~11× slower per call
  (207 s vs 19 s)** than the old 10 s budget — for no F1 gain. The old wall-clock
  cap was cutting the optimizer off early at essentially zero F1 cost.
- **Reproducibility goal: achieved** (output is now deterministic). But the
  chosen `MAX_OPT_PASSES=5` buys determinism at a large, unnecessary time cost.

**Decision: the set-level optimizer was removed entirely.** No-opt F1 (0.9276)
is within noise of every optimized variant (old wall-clock 0.9317, passes=5
0.9302, passes=2 0.9300), so the optimizer added no motif-similarity value while
adding runtime cost and ~185 lines of code. After removal, the benchmark is
byte-identical to the no-opt run (ALL F1 0.9276, 40/40 replicates), and the
pipeline ends deterministically at `palindromize → dedup`. (Re-run on all 500
to confirm at scale.)

### Correction to the earlier write-up

A previous version of this file reported the new commit **regressing F1 by
−0.21**. That was wrong — an artifact of an eval that (a) scored raw
`find_motifs` output without MicrobeMod's N-stripping/re-masking, (b) collapsed
multi-mod genomes to one nondeterministically-chosen mod-type, capping recall.
With the production-faithful harness there is **no F1 regression**.

_Data: `eval_output/results_v2/{A_noopt,B_oldopt_10s,C_newopt_5pass}/`
(each: `eval_summary.json`, `eval_results.csv`, `motifs/<rep>.json`)._
