#!/usr/bin/env python3
"""Check rolling/SIMD output equality and report repeated kernel timings."""

import argparse
import pathlib
import re
import statistics
import subprocess


def run(binary: pathlib.Path, reference: pathlib.Path, reads: pathlib.Path, mode: str):
    process = subprocess.run(
        [
            str(binary),
            "--reference",
            str(reference),
            "--reads",
            str(reads),
            "--mode",
            mode,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    match = re.search(r"seconds=([0-9.]+).*GCUPS=([0-9.]+)", process.stderr)
    if not match:
        raise RuntimeError(f"Could not parse timing output: {process.stderr.strip()}")
    return process.stdout, float(match.group(1)), float(match.group(2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--reference", required=True, type=pathlib.Path)
    parser.add_argument("--reads", required=True, type=pathlib.Path)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()

    if args.repeats <= 0:
        parser.error("repeats must be positive")
    binary = args.binary.resolve()
    reference = args.reference.resolve()
    reads = args.reads.resolve()

    rolling_csv, _, _ = run(binary, reference, reads, "rolling")
    simd_csv, _, _ = run(binary, reference, reads, "simd")
    if rolling_csv != simd_csv:
        raise SystemExit("ERROR: rolling and SIMD scores or coordinates differ")
    print("correctness=exact (rolling and SIMD CSV outputs match)")

    # Warm up both paths before recording measurements.
    run(binary, reference, reads, "rolling")
    run(binary, reference, reads, "simd")

    times = {"rolling": [], "simd": []}
    gcups = {"rolling": [], "simd": []}
    for repeat in range(args.repeats):
        modes = ("rolling", "simd") if repeat % 2 == 0 else ("simd", "rolling")
        for mode in modes:
            _, seconds, rate = run(binary, reference, reads, mode)
            times[mode].append(seconds)
            gcups[mode].append(rate)

    medians = {mode: statistics.median(values) for mode, values in times.items()}
    for mode in ("rolling", "simd"):
        print(
            f"{mode}: median_seconds={medians[mode]:.6f} "
            f"median_GCUPS={statistics.median(gcups[mode]):.3f} "
            f"runs=" + ",".join(f"{value:.6f}" for value in times[mode])
        )
    print(f"speedup_rolling_over_simd={medians['rolling'] / medians['simd']:.3f}x")


if __name__ == "__main__":
    main()
