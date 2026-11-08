# Benchmark report: Smith–Waterman implementations

Generated: 2026-11-08T15:57:26+01:00

## Machine and build

- Machine: Apple M5 Pro
- OS: macOS-26.7.2-arm64-arm-64bit
- C++ compiler: Apple clang version 21.0.0 (clang-2100.3.34.2)
- CMake: cmake version 4.4.5
- Build: CMake Release (`-O3 -DNDEBUG`)
- External implementation: [giorgapost/smith-waterman-algorithm](https://github.com/giorgapost/smith-waterman-algorithm), commit `98660d57f818a1d18a6129c6e66947d598acf878`, GPL-3.0
- Timed repetitions: 84 total runs (7 per profile and implementation); each mode is warmed up once before recording
- DP cell updates across timed runs: 44,454,144,000

## Method

The comparison uses the same generated reference/read files, order, and scoring (`match=2`, `mismatch=-1`, `gap=-2`). Each profile runs the external sequential full-matrix implementation and local `full`, `rolling`, and ARM NEON `simd` implementations. The harness checks all scores and aborts unless every local reference interval is one of the external implementation's optimal intervals.

`DP time` is measured around the cell computation loops. Local `dp_seconds` excludes scratch-buffer allocation; the external `calcCellsTime` excludes allocation of its DP matrix. The external program also traces every optimal alignment, so its total runtime is recorded separately and is not compared to local total time. The local `total_seconds` includes local scratch allocations but not FASTA input or CSV output.

The synthetic generator selects read substrings from a seeded random DNA reference and applies the specified number of substitutions. Variable-length profiles test SIMD lane masking. GCUPS is calculated as `(reference length × sum of read lengths) / DP seconds / 10^9`.

## Dataset profiles

| Profile | Reference | Reads | Read lengths | Substitutions/read | Seed | DP cells/run |
|---|---:|---:|---:|---:|---:|---:|
| small | 4,000 bp | 64 | 120–120 bp | 4 | 4101 | 30,720,000 |
| medium | 16,000 bp | 128 | 101–220 bp | 5 | 16021 | 329,312,000 |
| large | 32,000 bp | 192 | 120–280 bp | 7 | 32021 | 1,227,616,000 |

## Results

### Small profile

| Implementation | Median DP (s) | Std dev (s) | Min–max (s) | GCUPS | Speedup vs online sequential |
|---|---:|---:|---:|---:|---:|
| `online` | 0.098456 | 0.000407 | 0.098098–0.099259 | 0.312 | 1.000× |
| `full` | 0.080130 | 0.001014 | 0.079926–0.082283 | 0.383 | 1.229× |
| `rolling` | 0.098150 | 0.004285 | 0.091525–0.102115 | 0.313 | 1.003× |
| `simd` | 0.037319 | 0.000071 | 0.037157–0.037372 | 0.823 | 2.638× |

DP cells per run: 30,720,000. External sequential total time median (including traceback): 0.101473 s.

### Medium profile

| Implementation | Median DP (s) | Std dev (s) | Min–max (s) | GCUPS | Speedup vs online sequential |
|---|---:|---:|---:|---:|---:|
| `online` | 1.053870 | 0.003936 | 1.052900–1.063990 | 0.312 | 1.000× |
| `full` | 0.846144 | 0.001094 | 0.843731–0.846696 | 0.389 | 1.245× |
| `rolling` | 0.997020 | 0.049979 | 0.952373–1.084015 | 0.330 | 1.057× |
| `simd` | 0.504184 | 0.001911 | 0.501962–0.507522 | 0.653 | 2.090× |

DP cells per run: 329,312,000. External sequential total time median (including traceback): 1.074970 s.

### Large profile

| Implementation | Median DP (s) | Std dev (s) | Min–max (s) | GCUPS | Speedup vs online sequential |
|---|---:|---:|---:|---:|---:|
| `online` | 4.053480 | 0.310065 | 3.913450–4.788640 | 0.303 | 1.000× |
| `full` | 3.254567 | 0.090520 | 3.145944–3.381023 | 0.377 | 1.245× |
| `rolling` | 3.665711 | 0.105482 | 3.571416–3.864837 | 0.335 | 1.106× |
| `simd` | 2.165430 | 0.022346 | 2.146703–2.203435 | 0.567 | 1.872× |

DP cells per run: 1,227,616,000. External sequential total time median (including traceback): 4.137320 s.

## Observations

NEON SIMD had the highest measured DP throughput in each profile. Its speedup versus the external sequential implementation for small, medium, and large was 2.64×, 2.09×, 1.87×. The speedup narrows on larger and more variable-length input; no hardware performance counters were collected, so cache/bandwidth effects are plausible explanations but are not measured here.

The `full` implementation had lower DP time than `rolling` in all three profiles, even though its peak memory was much higher. The rolling-row algorithm is a memory optimization; it does not guarantee lower runtime on this CPU and input mix.

## Raw data and reproduction

Every measured run is in [`benchmark-raw-results.csv`](benchmark-raw-results.csv); the per-profile medians and spread are in [`benchmark-summary.csv`](benchmark-summary.csv). FASTA files are generated deterministically by `scripts/generate_benchmark.py` and can be recreated by `scripts/benchmark_suite.py`.

For interpretation and caveats, see the main [README](../README.md) and the separate learning notes at the workspace root.

## Peak resident memory

The following peak RSS measurements are one process run per profile/mode, collected with macOS `/usr/bin/time -l`. They include FASTA parsing and result storage. This table complements the algorithmic memory complexity; it is not part of the repeated timing sample.

| Profile | Online sequential | `full` | `rolling` | `simd` |
|---|---:|---:|---:|---:|
| small | 4.50 MiB | 12.72 MiB | 1.81 MiB | 2.02 MiB |
| medium | 21.97 MiB | 323.08 MiB | 2.53 MiB | 3.28 MiB |
| large | 53.11 MiB | 1283.91 MiB | 3.44 MiB | 4.97 MiB |

On the large profile, full-matrix RSS reached 1283.91 MiB, compared with 3.44 MiB for rolling rows and 4.97 MiB for SIMD. The full implementation stores score and two start coordinates in every `Cell`; this is substantially larger than a matrix of integer scores alone. RSS also reflects parsed inputs and allocator behavior.

Raw RSS samples: [`benchmark-memory.csv`](benchmark-memory.csv).
