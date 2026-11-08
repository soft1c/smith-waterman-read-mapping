#!/usr/bin/env python3
"""Compare local modes with the pinned public sequential implementation."""

import argparse
import csv
import pathlib
import re
import statistics
import subprocess


ONLINE_REPOSITORY = "https://github.com/giorgapost/smith-waterman-algorithm"
ONLINE_COMMIT = "98660d57f818a1d18a6129c6e66947d598acf878"

DRIVER_SOURCE = r'''#include "SequentialImplementation.h"
#include <iostream>

class BenchmarkSequential final : public SequentialImplementation {
public:
    using SequentialImplementation::SequentialImplementation;

    void printAlignmentSummary() const {
        std::cout << "ONLINE_ALIGNMENT_RESULTS_BEGIN\n";
        for (std::size_t i = 0; i < results.size(); ++i) {
            int best = 0;
            for (const auto& alignment : results[i]) {
                if (alignment.score > best) best = alignment.score;
            }
            std::cout << "SCORE," << i << ',' << best << '\n';
            for (const auto& alignment : results[i]) {
                if (alignment.score == best) {
                    std::cout << "ALIGN," << i << ',' << alignment.score << ','
                              << alignment.start << ',' << (alignment.stop + 1) << '\n';
                }
            }
        }
        std::cout << "ONLINE_ALIGNMENT_RESULTS_END\n";
    }
};

int main(int argc, char** argv) {
    try {
        BenchmarkSequential implementation(argc, argv);
        implementation.runAlgorithm();
        implementation.printStatistics();
        implementation.printAlignmentSummary();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
'''


def read_fasta(path: pathlib.Path):
    records = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            records.append([line[1:].strip(), ""])
        elif not records:
            raise ValueError(f"Sequence before FASTA header in {path}")
        else:
            records[-1][1] += "".join(line.split()).upper()
    if not records or any(not sequence for _, sequence in records):
        raise ValueError(f"Empty or missing FASTA sequence in {path}")
    return records


def compile_online(source: pathlib.Path, output: pathlib.Path, compiler: str):
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if commit != ONLINE_COMMIT:
        raise SystemExit(
            f"Online source must be checked out at {ONLINE_COMMIT}; found {commit}"
        )

    output.mkdir(parents=True, exist_ok=True)
    driver = output / "online_benchmark_driver.cpp"
    binary = output / "online_sw_sequential"
    driver.write_text(DRIVER_SOURCE)
    subprocess.run(
        [
            compiler,
            "-O3",
            "-DNDEBUG",
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-I",
            str(source / "src/include"),
            str(source / "src/classes/Framework.cpp"),
            str(source / "src/classes/SequentialImplementation.cpp"),
            str(driver),
            "-o",
            str(binary),
        ],
        check=True,
    )
    return binary


def run_local(binary: pathlib.Path, reference: pathlib.Path,
              reads: pathlib.Path, mode: str):
    process = subprocess.run(
        [str(binary), "--reference", str(reference), "--reads", str(reads), "--mode", mode],
        capture_output=True,
        text=True,
        check=True,
    )
    timings = re.search(r"seconds=([0-9.]+).*dp_seconds=([0-9.]+)", process.stderr)
    if not timings:
        raise RuntimeError(f"Could not parse local timing: {process.stderr.strip()}")
    rows = list(csv.DictReader(process.stdout.splitlines()))
    return rows, float(timings.group(1)), float(timings.group(2))


def run_online(binary: pathlib.Path, pair_file: pathlib.Path):
    process = subprocess.run(
        [
            str(binary),
            "-path",
            str(pair_file),
            "-id",
            "benchmark",
            "-match",
            "2",
            "-mismatch",
            "-1",
            "-gap",
            "-2",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    dp_time = re.search(
        r"E\) Total time of calculating cells: ([0-9.eE+-]+)", process.stdout
    )
    total_time = re.search(
        r"D\) Total time of program execution: ([0-9.eE+-]+)", process.stdout
    )
    if not dp_time or not total_time:
        raise RuntimeError(f"Could not parse online timing: {process.stdout[:1000]}")

    scores = {}
    intervals = {}
    for line in process.stdout.splitlines():
        if line.startswith("SCORE,"):
            _, index, score = line.split(",")
            scores[int(index)] = int(score)
        elif line.startswith("ALIGN,"):
            _, index, score, start, end = line.split(",")
            intervals.setdefault(int(index), set()).add(
                (int(score), int(start), int(end))
            )
    return scores, intervals, float(dp_time.group(1)), float(total_time.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--online-source", required=True, type=pathlib.Path)
    parser.add_argument("--reference", required=True, type=pathlib.Path)
    parser.add_argument("--reads", required=True, type=pathlib.Path)
    parser.add_argument("--work-dir", type=pathlib.Path, default=pathlib.Path("build/online-benchmark"))
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()

    if args.repeats <= 0:
        parser.error("repeats must be positive")
    binary = args.binary.resolve()
    source = args.online_source.resolve()
    reference_path = args.reference.resolve()
    reads_path = args.reads.resolve()
    work_dir = args.work_dir.resolve()
    online_binary = compile_online(source, work_dir, args.compiler)

    references = read_fasta(reference_path)
    if len(references) != 1:
        parser.error("the online project accepts a single reference sequence per read")
    reads = read_fasta(reads_path)
    reference = references[0][1]
    pair_file = work_dir / "online-input.txt"
    pair_file.write_text("".join(f"Q: {sequence}\nD: {reference}\n\n" for _, sequence in reads))

    online_scores, online_intervals, _, _ = run_online(online_binary, pair_file)
    local_results = {
        mode: run_local(binary, reference_path, reads_path, mode)[0]
        for mode in ("full", "rolling", "simd")
    }
    expected_scores = [online_scores[index] for index in range(len(reads))]
    for mode, rows in local_results.items():
        if [int(row["score"]) for row in rows] != expected_scores:
            raise SystemExit(f"ERROR: online and local {mode} scores differ")

    interval_matches = {
        mode: sum(
            (
                int(row["score"]),
                int(row["reference_start"]),
                int(row["reference_end"]),
            )
            in online_intervals.get(index, set())
            for index, row in enumerate(rows)
        )
        for mode, rows in local_results.items()
    }
    print(f"score_check=exact across online/full/rolling/simd ({len(reads)} reads)")
    for mode, count in interval_matches.items():
        print(f"{mode}_interval_in_online_optimal_set={count}/{len(reads)}")
        if count != len(reads):
            raise SystemExit(f"ERROR: {mode} intervals differ from the online optimal set")
    print(f"online_source={ONLINE_REPOSITORY}@{ONLINE_COMMIT}")

    # Warm up all versions before collecting the alternating measurements.
    run_online(online_binary, pair_file)
    for mode in ("full", "rolling", "simd"):
        run_local(binary, reference_path, reads_path, mode)

    cells = sum(len(sequence) * len(reference) for _, sequence in reads)
    timings = {
        key: []
        for key in ("online_dp", "online_total", "full_dp", "rolling_dp", "simd_dp")
    }
    raw_rows = []
    for repeat in range(args.repeats):
        sequence = ["online", "full", "rolling", "simd"]
        if repeat % 2:
            sequence.reverse()
        for mode in sequence:
            if mode == "online":
                _, _, dp_time, total_time = run_online(online_binary, pair_file)
                timings["online_dp"].append(dp_time)
                timings["online_total"].append(total_time)
                raw_rows.append({
                    "repeat": repeat + 1,
                    "mode": mode,
                    "dp_seconds": dp_time,
                    "total_seconds": total_time,
                    "gcups": cells / dp_time / 1e9,
                })
            else:
                _, total_time, dp_time = run_local(
                    binary, reference_path, reads_path, mode)
                timings[f"{mode}_dp"].append(dp_time)
                raw_rows.append({
                    "repeat": repeat + 1,
                    "mode": mode,
                    "dp_seconds": dp_time,
                    "total_seconds": total_time,
                    "gcups": cells / dp_time / 1e9,
                })

    raw_path = work_dir / "benchmark-results.csv"
    with raw_path.open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=("repeat", "mode", "dp_seconds", "total_seconds", "gcups"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(raw_rows)
    for row in raw_rows:
        print(
            f"run={row['repeat']:02d} mode={row['mode']} "
            f"dp_seconds={row['dp_seconds']:.6f} "
            f"total_seconds={row['total_seconds']:.6f} "
            f"GCUPS={row['gcups']:.3f}"
        )

    medians = {key: statistics.median(values) for key, values in timings.items()}
    for mode in ("online_dp", "full_dp", "rolling_dp", "simd_dp"):
        seconds = medians[mode]
        print(f"{mode}_median_seconds={seconds:.6f} GCUPS={cells / seconds / 1e9:.3f}")
    print(f"online_total_median_seconds={medians['online_total']:.6f}")
    for mode in ("full_dp", "rolling_dp", "simd_dp"):
        local_mode = mode.removesuffix("_dp")
        print(
            f"speedup_local_{local_mode}_vs_online="
            f"{medians['online_dp'] / medians[mode]:.3f}x"
        )
    print(f"dataset_dp_cells={cells} reads={len(reads)} reference_length={len(reference)}")
    print(f"raw_measurements_csv={raw_path}")


if __name__ == "__main__":
    main()
