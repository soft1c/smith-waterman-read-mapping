#!/usr/bin/env python3
"""Compare full, rolling, and compact traceback for one read."""

import argparse
import pathlib
import re
import statistics
import subprocess


def run(
    binary: pathlib.Path,
    reference: pathlib.Path,
    reads: pathlib.Path,
    mode: str,
    extra_options=(),
):
    process = subprocess.run(
        [
            str(binary),
            "--reference",
            str(reference),
            "--reads",
            str(reads),
            "--mode",
            mode,
            *extra_options,
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
        help="optional prior build; compared in --baseline-mode",
    )
    parser.add_argument(
        "--baseline-mode",
        choices=("trace", "wavefront"),
        default="trace",
        help="mode used for the optional baseline binary (default: trace)",
    )
    parser.add_argument(
        "--trace-only",
        action="store_true",
        help="compare only trace modes; useful when full-matrix mode would use too much memory",
    )
    parser.add_argument(
        "--include-wavefront",
        action="store_true",
        help="also benchmark and compare single-read tiled parallelism",
    )
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--tile-rows", type=int, default=0)
    parser.add_argument("--tile-columns", type=int, default=2048)
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("repeats must be positive")
    if args.trace_only and not args.baseline_binary:
        parser.error("--trace-only requires --baseline-binary")
    if args.baseline_mode == "wavefront" and not args.include_wavefront:
        parser.error("--baseline-mode wavefront requires --include-wavefront")

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
    if args.include_wavefront:
        modes = (*modes, "wavefront")
    wavefront_options = (
        "--threads",
        str(args.threads),
        "--tile-rows",
        str(args.tile_rows),
        "--tile-columns",
        str(args.tile_columns),
    )

    def run_mode(measured_binary, mode):
        return run(
            measured_binary,
            reference,
            reads,
            mode,
            wavefront_options if mode == "wavefront" else (),
        )

    outputs = {mode: run_mode(binary, mode)[0] for mode in modes}
    reference_mode = "trace" if args.trace_only else "rolling"
    if any(outputs[mode] != outputs[reference_mode] for mode in modes):
        raise SystemExit(f"ERROR: candidate outputs differ from {reference_mode}")
    if baseline_binary is not None:
        baseline_output = run_mode(baseline_binary, args.baseline_mode)[0]
        if baseline_output != outputs[reference_mode]:
            raise SystemExit(
                f"ERROR: baseline {args.baseline_mode} and candidate "
                f"{reference_mode} outputs differ"
            )
    if args.trace_only and args.include_wavefront:
        print("correctness=exact (trace and wavefront CSV outputs match)")
    elif args.trace_only:
        print("correctness=exact (candidate and baseline trace CSV outputs match)")
    elif args.include_wavefront:
        print("correctness=exact (full, rolling, trace, and wavefront CSV outputs match)")
    else:
        print("correctness=exact (full, rolling, and trace CSV outputs match)")

    for mode in modes:
        run_mode(binary, mode)
    measurements = [(mode, binary, mode) for mode in modes]
    baseline_label = f"baseline-{args.baseline_mode}"
    if baseline_binary is not None:
        run_mode(baseline_binary, args.baseline_mode)
        measurements.append((baseline_label, baseline_binary, args.baseline_mode))
    times = {label: [] for label, _, _ in measurements}
    gcups = {label: [] for label, _, _ in measurements}
    for repeat in range(args.repeats):
        offset = repeat % len(measurements)
        run_order = measurements[offset:] + measurements[:offset]
        for label, measured_binary, mode in run_order:
            _, seconds, rate = run_mode(measured_binary, mode)
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
    if args.include_wavefront:
        print(
            "speedup_trace_over_wavefront="
            f"{medians['trace'] / medians['wavefront']:.3f}x"
        )
    if baseline_binary is not None:
        print(
            f"speedup_candidate_{args.baseline_mode}_over_baseline_"
            f"{args.baseline_mode}="
            f"{medians[baseline_label] / medians[args.baseline_mode]:.3f}x"
        )


if __name__ == "__main__":
    main()
