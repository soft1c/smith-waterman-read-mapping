#!/usr/bin/env python3
"""Check rolling/trace/SIMD output equality and report repeated timings."""

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

    modes = ("rolling", "trace", "simd")
    outputs = {mode: run(binary, reference, reads, mode)[0] for mode in modes}
    if any(outputs[mode] != outputs["rolling"] for mode in modes[1:]):
        raise SystemExit("ERROR: rolling, trace, and SIMD CSV outputs differ")
    print("correctness=exact (rolling, trace, and SIMD CSV outputs match)")

    # Warm up each path before recording measurements.
    for mode in modes:
        run(binary, reference, reads, mode)

    times = {mode: [] for mode in modes}
    gcups = {mode: [] for mode in modes}
    for repeat in range(args.repeats):
        offset = repeat % len(modes)
        run_order = modes[offset:] + modes[:offset]
        for mode in run_order:
            _, seconds, rate = run(binary, reference, reads, mode)
            times[mode].append(seconds)
            gcups[mode].append(rate)

    medians = {mode: statistics.median(values) for mode, values in times.items()}
    for mode in modes:
        print(
            f"{mode}: median_seconds={medians[mode]:.6f} "
            f"median_GCUPS={statistics.median(gcups[mode]):.3f} "
            f"runs=" + ",".join(f"{value:.6f}" for value in times[mode])
        )
    for mode in modes[1:]:
        print(f"speedup_rolling_over_{mode}={medians['rolling'] / medians[mode]:.3f}x")


if __name__ == "__main__":
    main()
