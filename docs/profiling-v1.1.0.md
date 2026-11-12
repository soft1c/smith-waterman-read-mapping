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

## Next profiler iteration: v1.3.0

To resolve source frames in the Apple sample, a profiling-only Release build used `-O3 -DNDEBUG -g -fno-inline`. This changes performance substantially, so only the call graph was used to choose the next change. In the v1.2.0 capture, samples appeared in both zero-initializing the direction vector and the per-cell `std::max(initializer_list)` call. The optimized v1.2.0 sample otherwise attributed its frames to `main` because the compiler inlined the kernel.

The v1.3.0 `trace` implementation packs four directions into a row-local byte and writes that byte once. It allocates the direction array without zero-filling it, since the DP pass assigns every packed byte. The score maximum is also written as direct comparisons. This removes repeated read-modify-write direction updates and the separate buffer initialization pass.

| Single-read input | v1.2.0 trace median | v1.3.0 trace median | Speedup | Repetitions |
|---|---:|---:|---:|---:|
| Reference 100,000 bp; read 160 bp | 0.029166 s | 0.027973 s | 1.043× | 7 |
| Reference 100,000 bp; read 10,000 bp | 1.562736 s | 1.465844 s | 1.066× | 9 |

CSV output matched exactly between versions. The second profile's direction matrix still uses about 240 MiB; v1.3.0 reduces the work to initialize and update it, not the amount of traceback data stored. A separate odd-width reference check (12,001 bp, four 161 bp reads) also matched across rolling, trace, and SIMD.

The `-fno-inline` sample for v1.3.0 no longer showed the direction-vector zero-fill or initializer-list maximum as hot frames. Remaining source-mapped samples were concentrated around sequence-base access and score-row loads. Because this profile build is much slower than the optimized build, the 1.043× and 1.066× figures above come from separate optimized benchmark runs, not the profile capture.

Reproduce that source-mapped capture with:

```sh
cmake -S . -B build/profile-noinline -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g -fno-inline'
cmake --build build/profile-noinline -j
build/profile-noinline/sw_mapper --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --mode trace \
  > build/profiling/trace-v1.3-noinline.stdout.csv \
  2> build/profiling/trace-v1.3-noinline.stderr.txt &
sw_pid=$!
/usr/bin/sample "$sw_pid" 1 1 -mayDie -fullPaths \
  -file build/profiling/trace-v1.3-noinline.sample.txt
wait "$sw_pid"
```

Compare the v1.3.0 candidate against a v1.2.0 build with:

```sh
python3 scripts/benchmark_single.py --binary build/profile-release/sw_mapper \
  --baseline-binary build/releases/sw_mapper-v1.2.0 --trace-only \
  --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --repeats 9
```

## Next profiler iteration: v1.4.0

The v1.3.0 source-mapped sample concentrated on per-cell base access and score-row indexing. In v1.4.0, the read base and raw pointers to the reference and score rows are hoisted outside the cell loop. Direction bytes are assembled from four consecutive cells and written once per block, removing the shift counter and pack-flush branch from each cell. The resulting `sample` capture still maps the optimized kernel into `main` because the compiler inlines it; the `-fno-inline` source capture maps the hot samples to the per-cell lambda and recurrence. It does not provide hardware wait or cache counters.

| Single-read input | v1.3.0 trace median | v1.4.0 trace median | Speedup | Repetitions |
|---|---:|---:|---:|---:|
| Reference 100,000 bp; read 160 bp | 0.028194 s | 0.021950 s | 1.284× | 7 |
| Reference 100,000 bp; read 10,000 bp | 1.550714 s | 1.232614 s | 1.258× | 5 |

The benchmark compared outputs exactly. On the 160 bp case, v1.4.0 was 2.44× faster than rolling v1.1.0 and 3.07× faster than full v1.0.0. The odd-width reference check also still matched rolling and SIMD. This run confirms a single-read speedup on both tested read lengths; it does not establish the same gain for other sequence sizes or scoring parameters.

The optimized 1-second Apple `sample` recorded 857 stack samples, attributed to the inlined `main` call site. The run's reported wall time was 1.303 s and DP time 1.303 s. A source-mapped `-fno-inline` capture completed in 1.527 s and showed the recurrence as the remaining hot region. These samples still do not say whether the recurrence waits on memory or execution units.

Compare v1.4.0 with v1.3.0 using:

```sh
python3 scripts/benchmark_single.py --binary build/profile-release/sw_mapper \
  --baseline-binary build/releases/sw_mapper-v1.3.0 \
  --reference build/single-read-data/reference.fa \
  --reads build/single-read-data/reads.fa --repeats 7
```

This is a CPU stack sample, not a hardware-counter trace. It cannot tell whether the hot instructions are limited by cache misses, memory bandwidth, branch prediction, or instruction throughput. The active developer directory is `/Library/Developer/CommandLineTools`; `xcrun xctrace` reports that `xctrace` is unavailable, and no full Xcode/Instruments installation was found. An Instruments Time Profiler and memory/counter trace should follow when that toolchain is available. Until then, this report does not claim a measured memory-stall bottleneck.

## Single-read parallel iteration: v1.5.0

The v1.4.0 profile left the scalar Smith–Waterman recurrence as the main cost. A read's DP grid has dependencies from the cell above and the cell to the left, so individual cells cannot be split arbitrarily between threads. v1.5.0 divides the grid into tiles and computes tiles along anti-diagonals. A tile starts after its north and west neighbors finish; tiles on the same anti-diagonal run concurrently. The main thread also claims tiles from the shared queue, while a persistent worker pool handles the remaining work for that alignment.

Each tile retains only its bottom and right score borders. The packed 2-bit traceback directions are still recorded for each cell, and tile widths must be a multiple of four so neighboring workers never write the same packed byte. Tile-best endpoints are reduced by score and then earliest row/column, preserving the sequential `trace` tie behavior. Use `--mode wavefront`; `--threads 0` selects `std::thread::hardware_concurrency()`, and `--tile-rows 0` picks 8 rows for reads up to 512 bp or 64 rows for longer reads. `--tile-columns` defaults to 2048 and must be a multiple of four.

The optimized build used Apple clang 21.0.0 with `-O3 -DNDEBUG -g` on an Apple M5 Pro / ARM64 running macOS 26.7.2. Single-read measurements compared the tagged v1.4.0 `trace` binary with v1.5.0 `wavefront --threads 0`; runs were interleaved, outputs were checked for exact CSV equality, and the reported values are medians:

| Input | v1.4.0 `trace` | v1.5.0 `wavefront` | Speedup | Runs |
|---|---:|---:|---:|---:|
| Reference 100,000 bp; read 160 bp | 0.022558 s | 0.007018 s | 3.21× | 9 |
| Reference 100,000 bp; read 10,000 bp | 1.288040 s | 0.155313 s | 8.29× | 7 |

The 160 bp workload used seed `21021`; the 10,000 bp workload used seed `21022`. The speedup depends on the machine's core count and read/reference shape. In a separate eight-case randomized check, `wavefront` matched `trace` with varied odd sequence lengths, tile boundaries, and scoring parameters.

Apple `/usr/bin/sample` captured a sustained run of twelve 10,000 bp alignments, processed one after another (the workload repeats the single-read kernel; reads are not parallelized with each other). Of 732 main-thread samples, 551 were in the tile-computation lambda and 139 in the condition-variable wait; the main thread now performs DP work instead of spending most of the run waiting for workers. Worker samples showed the same tile kernel and idle waits between anti-diagonals. This points to scheduling/barrier synchronization as the next software target; the stack sampler still cannot identify memory stalls or cache behavior. `xcrun xctrace` remains unavailable in the installed Command Line Tools.

On the 10,000 bp read, `/usr/bin/time -l` reported 241.1 MiB maximum RSS for v1.4.0 and 303.5 MiB for v1.5.0 wavefront (about 62.4 MiB more). The extra memory stores tile borders alongside the packed traceback matrix. The tested read completed about 8.3× faster with that memory increase.

Run the candidate with:

```sh
cmake -S . -B build/profile-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g'
cmake --build build/profile-release -j
build/profile-release/sw_mapper --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --mode wavefront --threads 0
```

## Synchronization iteration: v1.6.0

The v1.5.0 sampler showed the coordinator waiting for the whole worker group after every anti-diagonal. v1.6.0 replaces those global barriers with per-tile dependency counts and a bounded ready queue. A tile is queued as soon as its north and west neighbors finish; workers and the main thread take ready tiles directly. This lets a tile start without waiting for unrelated work elsewhere on the previous anti-diagonal.

The queue stores tile IDs in a preallocated ring, and each tile is enqueued once. Atomic dependency counters publish score borders from both predecessors; the queue mutex and condition variable protect ready work and wake idle threads. The persistent worker pool and main-thread participation from v1.5.0 remain in place.

The optimized v1.6.0 build was measured on the same Apple M5 Pro / Apple clang 21.0.0 setup and single-read inputs as v1.5.0. Runs of v1.5.0 and the candidate were interleaved; all outputs matched the v1.4.0 `trace` CSV exactly:

| Input | v1.5.0 `wavefront` | v1.6.0 `wavefront` | Speedup | Runs |
|---|---:|---:|---:|---:|
| Reference 100,000 bp; read 160 bp | 0.007426 s | 0.003805 s | 1.95× | 9 |
| Reference 100,000 bp; read 10,000 bp | 0.148622 s | 0.118705 s | 1.25× | 7 |

An eight-case randomized check also matched `trace` across varied odd sequence lengths, tile boundaries, scoring parameters, and thread counts. Peak RSS on the 10,000 bp read was effectively unchanged: 303.5 MiB for v1.5.0 and 303.6 MiB for v1.6.0.

Apple `sample` captured the same sustained workload of twelve 10,000 bp alignments processed one after another. In the v1.6.0 main-thread profile, 31 of 668 samples were waiting on the condition variable (4.6%), compared with 139 of 732 (19.0%) in v1.5.0. Only a few samples showed queue-mutex contention. Most samples stayed in the ready-tile processing loop. This supports that the repeated global wait was reduced; the remaining CPU work is in tile processing and queue handling. The capture is still a stack sample and does not provide cache or memory-stall counters. The profile is saved locally as `build/profiling/wavefront-v1.6-sample.txt`.

Reproduce the v1.5.0-to-v1.6.0 benchmark with:

```sh
python3 scripts/benchmark_single.py --binary build/profile-release/sw_mapper \
  --baseline-binary build/releases/sw_mapper-v1.5.0 \
  --baseline-mode wavefront --trace-only --include-wavefront \
  --reference build/single-profile-data/reference.fa \
  --reads build/single-profile-data/reads.fa --repeats 7
```

## Reproduction

Use the same generated data, optimized build flags, and sampling interval for later versions. Record the complete binary version, compiler, macOS version, dataset parameters, wall time, DP time, and profile file. Compare each single-read optimization against rolling v1.1.0 for output equality and repeated runtime; keep the sampler capture separate from unprofiled timing runs because profiling adds overhead. In v1.6.0, tile recurrence and ready-queue overhead remain software targets; hardware counters are still needed to classify any memory stalls.
