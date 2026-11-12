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
    parser.add_argument(
        "--baseline-binary",
        type=pathlib.Path,
        help="optional prior build; its trace mode is compared with the candidate",
    )
    parser.add_argument(
        "--trace-only",
        action="store_true",
        help="compare only trace modes; useful when full-matrix mode would use too much memory",
    )
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("repeats must be positive")
    if args.trace_only and not args.baseline_binary:
        parser.error("--trace-only requires --baseline-binary")

    binary = args.binary.resolve()
    baseline_binary = (
        args.baseline_binary.resolve() if args.baseline_binary else None
    )
    reference = args.reference.resolve()
    reads = args.reads.resolve()
    read_count = sum(line.startswith(">") for line in reads.read_text().splitlines())
    if read_count != 1:
        parser.error(f"provide exactly one read; found {read_count}")

    modes = ("trace",) if args.trace_only else ("full", "rolling", "trace")
    outputs = {mode: run(binary, reference, reads, mode)[0] for mode in modes}
    if not args.trace_only and any(
        outputs[mode] != outputs["rolling"] for mode in ("full", "trace")
    ):
        raise SystemExit("ERROR: full, rolling, and trace outputs differ")
    if baseline_binary is not None:
        baseline_output = run(baseline_binary, reference, reads, "trace")[0]
        if baseline_output != outputs["trace"]:
            raise SystemExit("ERROR: prior and candidate trace outputs differ")
    print(
        "correctness=exact (candidate and baseline trace CSV outputs match)"
        if args.trace_only
        else "correctness=exact (full, rolling, and trace CSV outputs match)"
    )

    for mode in modes:
        run(binary, reference, reads, mode)
    measurements = [(mode, binary, mode) for mode in modes]
    if baseline_binary is not None:
        run(baseline_binary, reference, reads, "trace")
        measurements.append(("baseline-trace", baseline_binary, "trace"))
    times = {label: [] for label, _, _ in measurements}
    gcups = {label: [] for label, _, _ in measurements}
    for repeat in range(args.repeats):
        offset = repeat % len(measurements)
        run_order = measurements[offset:] + measurements[:offset]
        for label, measured_binary, mode in run_order:
            _, seconds, rate = run(measured_binary, reference, reads, mode)
            times[label].append(seconds)
            gcups[label].append(rate)

    medians = {mode: statistics.median(samples) for mode, samples in times.items()}
    for mode in times:
        print(
            f"{mode}: median_seconds={medians[mode]:.6f} "
            f"median_GCUPS={statistics.median(gcups[mode]):.3f} "
            f"runs=" + ",".join(f"{value:.6f}" for value in times[mode])
        )
    if not args.trace_only:
        for mode in ("rolling", "trace"):
            print(f"speedup_{mode}_over_full={medians['full'] / medians[mode]:.3f}x")
        print(f"speedup_rolling_over_trace={medians['rolling'] / medians['trace']:.3f}x")
    if baseline_binary is not None:
        print(
            "speedup_candidate_trace_over_baseline_trace="
            f"{medians['baseline-trace'] / medians['trace']:.3f}x"
        )


if __name__ == "__main__":
    main()
