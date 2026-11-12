# Apple CPU profile: rolling baseline v1.1.0

## Capture

- Machine: Apple M5 Pro, ARM64
- OS: macOS 26.7.2
- Compiler: Apple clang 21.0.0
- Build: Release, `-O3 -DNDEBUG -g` (debug symbols enabled for source lines)
- Dataset: deterministic random reference, 100,000 bp; 256 reads, 160 bp each; seed `21021`
- Work: 4.096 billion DP cell updates
- Capture: Apple's `/usr/bin/sample`, 8 seconds, 1 ms sampling interval
- Result file: `build/profiling/rolling-v1.1.0.sample.txt` (local build output, not tracked)

The capture command was:

```sh
cmake -S . -B build/profile-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g'
cmake --build build/profile-release -j
python3 scripts/generate_benchmark.py --output-dir build/profile-data \
  --reference-length 100000 --read-count 256 --read-length 160 \
  --substitutions 5 --seed 21021
build/profile-release/sw_mapper --reference build/profile-data/reference.fa \
  --reads build/profile-data/reads.fa --mode rolling \
  > build/profiling/rolling.stdout.csv 2> build/profiling/rolling.stderr.txt &
sw_pid=$!
/usr/bin/sample "$sw_pid" 8 1 \
  -file build/profiling/rolling-v1.1.0.sample.txt -fullPaths
wait "$sw_pid"
```

## Findings

The run took 12.884 s wall time, of which 12.878 s was reported inside the DP timer (99.96%). The 8-second profile recorded 6,846 stack samples. They all led through `main` into `smith_waterman_rolling`; there was no notable time in FASTA parsing, output, or allocation in this capture.

| Sample location | Samples | Share |
|---|---:|---:|
| `src/main.cpp:198` — cell selection and start-coordinate propagation | 5,442 | 79.5% |
| `src/main.cpp:204` — updating the best alignment | 715 | 10.4% |
| `src/main.cpp:200` — substitution score | 435 | 6.4% |
| `src/main.cpp:194` — inner-loop overhead and other mapped instructions | 229 | 3.3% |

The profile identified per-cell work in the rolling kernel, especially carrying `Cell` start coordinates through the recurrence. The single-read v1.2.0 iteration below removes that coordinate propagation from the score loop.

## Single-read iteration: v1.2.0

The new `trace` mode keeps only the score values for the previous and current row. It stores the winning direction for each cell in a packed 2-bit matrix, then walks back from the highest-scoring cell to recover the start coordinates. This avoids reading and writing two `size_t` start coordinates on every DP cell. Direction storage takes about one quarter of a byte per DP cell, so it trades memory for less work in the hot loop.

The benchmark used one 160 bp read against a 100,000 bp reference, with five substitutions and seed `21021`. The optimized build used `-O3 -DNDEBUG -g`; the script checked exact CSV equality, warmed all three modes once, then ran seven measurements in rotated order.

| Mode | Median wall time | GCUPS | Timed run range |
|---|---:|---:|---|
| `full` v1.0.0 | 0.063017 s | 0.254 | 0.061450–0.064454 s |
| `rolling` v1.1.0 | 0.051301 s | 0.312 | 0.050085–0.051780 s |
| `trace` v1.2.0 | 0.028445 s | 0.562 | 0.026809–0.029458 s |

On this input, rolling was **1.23× faster** than full; trace was **2.22× faster** than full and **1.80× faster** than rolling. All three CSV outputs matched exactly. This is one fixed workload, so it does not establish the same speedup for longer reads. The packed direction matrix grows with both read and reference length: profiling one 10,000 bp read against a 100,000 bp reference recorded 240.9 MiB peak physical footprint.

For the 160 bp read, one `/usr/bin/time -l` process run reported 6.64 MiB peak RSS for rolling and 6.55 MiB for `trace`. This particular shape makes the packed direction data small enough that its peak memory was close to, and slightly below, the two rolling `Cell` rows. Longer reads increase `trace` memory linearly with the read length.

The 1-second Apple `sample` capture for that longer read recorded 856 stack samples, but the optimized build attributed them to the `main` call site rather than resolving the hot loop to source lines. It confirms the run stayed CPU-active; it does not identify cache misses or memory stalls. Hardware-counter analysis still requires Instruments, which is unavailable in the current Command Line Tools installation.

Capture the v1.2.0 sample with:

```sh
python3 scripts/generate_benchmark.py --output-dir build/single-profile-data \
  --reference-length 100000 --read-count 1 --read-length 10000 \
  --substitutions 12 --seed 21022
build/profile-release/sw_mapper --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --mode trace \
  > build/profiling/trace-v1.2.stdout.csv 2> build/profiling/trace-v1.2.stderr.txt &
sw_pid=$!
/usr/bin/sample "$sw_pid" 1 1 \
  -file build/profiling/trace-v1.2.sample.txt -fullPaths
wait "$sw_pid"
```

Reproduce the single-read benchmark with:

```sh
python3 scripts/generate_benchmark.py --output-dir build/single-read-data \
  --reference-length 100000 --read-count 1 --read-length 160 \
  --substitutions 5 --seed 21021
python3 scripts/benchmark_single.py --binary build/profile-release/sw_mapper \
  --reference build/single-read-data/reference.fa \
  --reads build/single-read-data/reads.fa --repeats 7
```

This is a CPU stack sample, not a hardware-counter trace. It cannot tell whether the hot instructions are limited by cache misses, memory bandwidth, branch prediction, or instruction throughput. The active developer directory is `/Library/Developer/CommandLineTools`; `xcrun xctrace` reports that `xctrace` is unavailable, and no full Xcode/Instruments installation was found. An Instruments Time Profiler and memory/counter trace should follow when that toolchain is available. Until then, this report does not claim a measured memory-stall bottleneck.

## Reproduction

Use the same generated data, optimized build flags, and sampling interval for later versions. Record the complete binary version, compiler, macOS version, dataset parameters, wall time, DP time, and profile file. Compare each single-read optimization against rolling v1.1.0 for output equality and repeated runtime; keep the sampler capture separate from unprofiled timing runs because profiling adds overhead. The next profiler pass should use Instruments once available to inspect the single-read kernel and hardware memory counters.
