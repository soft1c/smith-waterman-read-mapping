#!/usr/bin/env python3
"""Compare full, rolling, and compact traceback for one read."""

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
    read_count = sum(line.startswith(">") for line in reads.read_text().splitlines())
    if read_count != 1:
        parser.error(f"provide exactly one read; found {read_count}")

    modes = ("full", "rolling", "trace")
    outputs = {mode: run(binary, reference, reads, mode)[0] for mode in modes}
    if any(outputs[mode] != outputs["rolling"] for mode in ("full", "trace")):
        raise SystemExit("ERROR: full, rolling, and trace outputs differ")
    print("correctness=exact (full, rolling, and trace CSV outputs match)")

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

    medians = {mode: statistics.median(samples) for mode, samples in times.items()}
    for mode in modes:
        print(
            f"{mode}: median_seconds={medians[mode]:.6f} "
            f"median_GCUPS={statistics.median(gcups[mode]):.3f} "
            f"runs=" + ",".join(f"{value:.6f}" for value in times[mode])
        )
    for mode in ("rolling", "trace"):
        print(f"speedup_{mode}_over_full={medians['full'] / medians[mode]:.3f}x")
    print(f"speedup_rolling_over_trace={medians['rolling'] / medians['trace']:.3f}x")


if __name__ == "__main__":
    main()
