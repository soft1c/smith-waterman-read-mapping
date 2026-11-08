#!/usr/bin/env python3
"""Run repeated large-data comparisons and write a report plus raw CSV."""

import argparse
import csv
import datetime
import json
import pathlib
import platform
import statistics
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
EXTERNAL_COMMIT = "98660d57f818a1d18a6129c6e66947d598acf878"
PROFILES = (
    {
        "name": "small",
        "reference_length": 4_000,
        "read_count": 64,
        "min_read_length": 120,
        "max_read_length": 120,
        "substitutions": 4,
        "seed": 4101,
    },
    {
        "name": "medium",
        "reference_length": 16_000,
        "read_count": 128,
        "min_read_length": 100,
        "max_read_length": 220,
        "substitutions": 5,
        "seed": 16021,
    },
    {
        "name": "large",
        "reference_length": 32_000,
        "read_count": 192,
        "min_read_length": 120,
        "max_read_length": 280,
        "substitutions": 7,
        "seed": 32021,
    },
)


def command_output(command):
    return subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout.strip()


def fasta_lengths(path: pathlib.Path):
    lengths = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            lengths.append(0)
        elif line.strip():
            lengths[-1] += len("".join(line.split()))
    return lengths


def write_summary(profile_results, report_path, raw_path, metadata, repeats):
    raw_fields = (
        "profile",
        "reference_length",
        "read_count",
        "min_read_length",
        "max_read_length",
        "mean_read_length",
        "dp_cells",
        "repeat",
        "mode",
        "dp_seconds",
        "total_seconds",
        "gcups",
    )
    raw_rows = []
    summary_rows = []
    total_timed_cells = 0
    for profile, _dataset_dir, lengths, csv_path in profile_results:
        with csv_path.open(newline="") as input_file:
            rows = list(csv.DictReader(input_file))
        cell_count = sum(length * profile["reference_length"] for length in lengths)
        total_timed_cells += cell_count * 4 * repeats
        mean_length = statistics.mean(lengths)
        for row in rows:
            raw_rows.append({
                "profile": profile["name"],
                "reference_length": profile["reference_length"],
                "read_count": len(lengths),
                "min_read_length": min(lengths),
                "max_read_length": max(lengths),
                "mean_read_length": f"{mean_length:.2f}",
                "dp_cells": cell_count,
                "repeat": row["repeat"],
                "mode": row["mode"],
                "dp_seconds": row["dp_seconds"],
                "total_seconds": row["total_seconds"],
                "gcups": row["gcups"],
            })

        for mode in ("online", "full", "rolling", "simd"):
            selected = [row for row in rows if row["mode"] == mode]
            dp_times = [float(row["dp_seconds"]) for row in selected]
            total_times = [float(row["total_seconds"]) for row in selected]
            median_dp = statistics.median(dp_times)
            summary_rows.append({
                "profile": profile["name"],
                "mode": mode,
                "median_dp_seconds": median_dp,
                "stdev_dp_seconds": statistics.stdev(dp_times) if len(dp_times) > 1 else 0.0,
                "min_dp_seconds": min(dp_times),
                "max_dp_seconds": max(dp_times),
                "median_total_seconds": statistics.median(total_times),
                "median_gcups": cell_count / median_dp / 1e9,
                "speedup_over_online": (
                    statistics.median(
                        float(row["dp_seconds"])
                        for row in rows if row["mode"] == "online"
                    ) / median_dp
                ),
            })

    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=raw_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(raw_rows)

    summary_path = raw_path.with_name("benchmark-summary.csv")
    with summary_path.open("w", newline="") as output:
        fields = tuple(summary_rows[0].keys())
        writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "# Benchmark report: Smith–Waterman implementations",
        "",
        f"Generated: {datetime.datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Machine and build",
        "",
        f"- Machine: {metadata['machine']}",
        f"- OS: {metadata['platform']}",
        f"- C++ compiler: {metadata['compiler']}",
        f"- CMake: {metadata['cmake']}",
        "- Build: CMake Release (`-O3 -DNDEBUG`)",
        f"- External implementation: [giorgapost/smith-waterman-algorithm]({metadata['external_url']}), commit `{EXTERNAL_COMMIT}`, GPL-3.0",
        f"- Timed repetitions: {len(profile_results) * 4 * repeats} total runs ({repeats} per profile and implementation); each mode is warmed up once before recording",
        f"- DP cell updates across timed runs: {total_timed_cells:,}",
        "",
        "## Method",
        "",
        "The comparison uses the same generated reference/read files, order, and scoring (`match=2`, `mismatch=-1`, `gap=-2`). Each profile runs the external sequential full-matrix implementation and local `full`, `rolling`, and ARM NEON `simd` implementations. The harness checks all scores and aborts unless every local reference interval is one of the external implementation's optimal intervals.",
        "",
        "`DP time` is measured around the cell computation loops. Local `dp_seconds` excludes scratch-buffer allocation; the external `calcCellsTime` excludes allocation of its DP matrix. The external program also traces every optimal alignment, so its total runtime is recorded separately and is not compared to local total time. The local `total_seconds` includes local scratch allocations but not FASTA input or CSV output.",
        "",
        "The synthetic generator selects read substrings from a seeded random DNA reference and applies the specified number of substitutions. Variable-length profiles test SIMD lane masking. GCUPS is calculated as `(reference length × sum of read lengths) / DP seconds / 10^9`.",
        "",
        "## Dataset profiles",
        "",
        "| Profile | Reference | Reads | Read lengths | Substitutions/read | Seed | DP cells/run |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for profile, _, lengths, _ in profile_results:
        cell_count = profile["reference_length"] * sum(lengths)
        lines.append(
            f"| {profile['name']} | {profile['reference_length']:,} bp | "
            f"{len(lengths)} | {min(lengths)}–{max(lengths)} bp | "
            f"{profile['substitutions']} | {profile['seed']} | {cell_count:,} |"
        )

    lines += ["", "## Results", ""]
    for profile, _, lengths, _ in profile_results:
        lines += [
            f"### {profile['name'].capitalize()} profile",
            "",
            "| Implementation | Median DP (s) | Std dev (s) | Min–max (s) | GCUPS | Speedup vs online sequential |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        cell_count = profile["reference_length"] * sum(lengths)
        for row in summary_rows:
            if row["profile"] != profile["name"]:
                continue
            lines.append(
                f"| `{row['mode']}` | {row['median_dp_seconds']:.6f} | "
                f"{row['stdev_dp_seconds']:.6f} | "
                f"{row['min_dp_seconds']:.6f}–{row['max_dp_seconds']:.6f} | "
                f"{row['median_gcups']:.3f} | {row['speedup_over_online']:.3f}× |"
            )
        online_total = statistics.median(
            float(row["total_seconds"])
            for row in raw_rows
            if row["profile"] == profile["name"] and row["mode"] == "online"
        )
        lines += [
            "",
            f"DP cells per run: {cell_count:,}. External sequential total time median (including traceback): {online_total:.6f} s.",
            "",
        ]

    simd_speedups = [
        row["speedup_over_online"]
        for row in summary_rows
        if row["mode"] == "simd"
    ]
    lines += [
        "## Observations",
        "",
        "NEON SIMD had the highest measured DP throughput in each profile. Its speedup versus the external sequential implementation for small, medium, and large was "
        + ", ".join(f"{value:.2f}×" for value in simd_speedups)
        + ". The speedup narrows on larger and more variable-length input; no hardware performance counters were collected, so cache/bandwidth effects are plausible explanations but are not measured here.",
        "",
        "The `full` implementation had lower DP time than `rolling` in all three profiles, even though its peak memory was much higher. The rolling-row algorithm is a memory optimization; it does not guarantee lower runtime on this CPU and input mix.",
        "",
    ]

    lines += [
        "## Raw data and reproduction",
        "",
        f"Every measured run is in [`{raw_path.name}`]({raw_path.name}); the per-profile medians and spread are in [`{summary_path.name}`]({summary_path.name}). FASTA files are generated deterministically by `scripts/generate_benchmark.py` and can be recreated by `scripts/benchmark_suite.py`.",
        "",
        "For interpretation and caveats, see the main [README](../README.md) and the separate learning notes at the workspace root.",
        "",
    ]
    memory_path = raw_path.with_name("benchmark-memory.csv")
    if memory_path.exists():
        with memory_path.open(newline="") as memory_file:
            memory_rows = list(csv.DictReader(memory_file))
        lines += [
            "## Peak resident memory",
            "",
            "The following peak RSS measurements are one process run per profile/mode, collected with macOS `/usr/bin/time -l`. They include FASTA parsing and result storage. This table complements the algorithmic memory complexity; it is not part of the repeated timing sample.",
            "",
            "| Profile | Online sequential | `full` | `rolling` | `simd` |",
            "|---|---:|---:|---:|---:|",
        ]
        for profile, *_ in profile_results:
            values = {
                row["mode"]: float(row["peak_rss_mib"])
                for row in memory_rows
                if row["profile"] == profile["name"]
            }
            lines.append(
                f"| {profile['name']} | {values['online']:.2f} MiB | "
                f"{values['full']:.2f} MiB | {values['rolling']:.2f} MiB | "
                f"{values['simd']:.2f} MiB |"
            )
        large_memory = {
            row["mode"]: float(row["peak_rss_mib"])
            for row in memory_rows
            if row["profile"] == "large"
        }
        lines += [
            "",
            f"On the large profile, full-matrix RSS reached {large_memory['full']:.2f} MiB, compared with {large_memory['rolling']:.2f} MiB for rolling rows and {large_memory['simd']:.2f} MiB for SIMD. The full implementation stores score and two start coordinates in every `Cell`; this is substantially larger than a matrix of integer scores alone. RSS also reflects parsed inputs and allocator behavior.",
            "",
            f"Raw RSS samples: [`{memory_path.name}`]({memory_path.name}).",
            "",
        ]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines))
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--online-source", required=True, type=pathlib.Path)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--output-dir", type=pathlib.Path, default=ROOT / "build/benchmark-suite")
    parser.add_argument("--report", type=pathlib.Path, default=ROOT / "docs/benchmark-results.md")
    parser.add_argument("--raw-csv", type=pathlib.Path, default=ROOT / "docs/benchmark-raw-results.csv")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument(
        "--summarize-existing",
        action="store_true",
        help="create report from existing profile data and measurements without rerunning them",
    )
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("repeats must be positive")

    output_dir = args.output_dir.resolve()
    source = args.online_source.resolve()
    report_path = args.report.resolve()
    raw_path = args.raw_csv.resolve()
    binary = ROOT / "build/sw_mapper"
    if not args.skip_build:
        subprocess.run(
            ["cmake", "-S", str(ROOT), "-B", str(ROOT / "build"), "-DCMAKE_BUILD_TYPE=Release"],
            check=True,
        )
        subprocess.run(["cmake", "--build", str(ROOT / "build"), "--parallel"], check=True)

    metadata = {
        "machine": command_output(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "platform": platform.platform(),
        "compiler": command_output(["clang++", "--version"]).splitlines()[0],
        "cmake": command_output(["cmake", "--version"]).splitlines()[0],
        "external_url": "https://github.com/giorgapost/smith-waterman-algorithm",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_results = []
    if args.summarize_existing:
        for profile in PROFILES:
            profile_dir = output_dir / profile["name"]
            dataset_dir = profile_dir / "data"
            csv_path = profile_dir / "online-build/benchmark-results.csv"
            if not (dataset_dir / "reads.fa").exists() or not csv_path.exists():
                raise SystemExit(f"Missing existing dataset or measurements for {profile['name']}")
            profile_results.append(
                (profile, dataset_dir, fasta_lengths(dataset_dir / "reads.fa"), csv_path)
            )
        summary_path = write_summary(
            profile_results, report_path, raw_path, metadata, args.repeats
        )
        metadata_path = output_dir / "environment.json"
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"report={report_path}")
        print(f"raw_csv={raw_path}")
        print(f"summary_csv={summary_path}")
        print(f"environment_json={metadata_path}")
        return

    for profile in PROFILES:
        profile_dir = output_dir / profile["name"]
        dataset_dir = profile_dir / "data"
        generator = [
            sys.executable,
            str(ROOT / "scripts/generate_benchmark.py"),
            "--output-dir",
            str(dataset_dir),
            "--reference-length",
            str(profile["reference_length"]),
            "--read-count",
            str(profile["read_count"]),
            "--min-read-length",
            str(profile["min_read_length"]),
            "--max-read-length",
            str(profile["max_read_length"]),
            "--substitutions",
            str(profile["substitutions"]),
            "--seed",
            str(profile["seed"]),
        ]
        subprocess.run(generator, check=True, cwd=ROOT)
        benchmark = [
            sys.executable,
            str(ROOT / "scripts/benchmark_online.py"),
            "--binary",
            str(binary),
            "--online-source",
            str(source),
            "--reference",
            str(dataset_dir / "reference.fa"),
            "--reads",
            str(dataset_dir / "reads.fa"),
            "--work-dir",
            str(profile_dir / "online-build"),
            "--repeats",
            str(args.repeats),
        ]
        subprocess.run(benchmark, check=True, cwd=ROOT)
        lengths = fasta_lengths(dataset_dir / "reads.fa")
        profile_results.append(
            (profile, dataset_dir, lengths, profile_dir / "online-build/benchmark-results.csv")
        )

    memory_path = raw_path.with_name("benchmark-memory.csv")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/measure_memory.py"),
            "--binary",
            str(binary),
            "--suite-dir",
            str(output_dir),
            "--output",
            str(memory_path),
        ],
        check=True,
        cwd=ROOT,
    )
    summary_path = write_summary(profile_results, report_path, raw_path, metadata, args.repeats)
    metadata_path = output_dir / "environment.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"report={report_path}")
    print(f"raw_csv={raw_path}")
    print(f"summary_csv={summary_path}")
    print(f"environment_json={metadata_path}")


if __name__ == "__main__":
    main()
