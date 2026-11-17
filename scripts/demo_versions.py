#!/usr/bin/env python3
"""Build and run every recorded Smith-Waterman iteration on the same FASTA input.

Historical versions are compiled from their Git commits without checking them
out. The v1.7 batch candidate and final auto dispatcher are compiled from the
current working-tree source because those iterations have no Git tags.
"""

import argparse
import csv
import datetime
import pathlib
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import time


ROOT = pathlib.Path(__file__).resolve().parents[1]
VERSIONS = (
    ("v1.0.0", "bc38b86", "full"),
    ("v1.1.0", "bc38b86", "rolling"),
    ("v1.2.0", "v1.2.0", "trace"),
    ("v1.3.0", "v1.3.0", "trace"),
    ("v1.4.0", "v1.4.0", "trace"),
    ("v1.5.0", "v1.5.0", "wavefront"),
    ("v1.6.0", "v1.6.0", "wavefront"),
    ("v1.7.0-batch", "working-tree", "batch"),
    ("final-auto", "working-tree", "auto"),
)
BUILD_FLAGS = ("-O3", "-DNDEBUG", "-g", "-std=c++17", "-pthread")


def run_text(command, cwd=ROOT):
    completed = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, check=True
    )
    return completed.stdout.strip()


def fasta_sequences(path):
    sequences = []
    current = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        elif line.strip():
            current.append("".join(line.split()))
    if current:
        sequences.append("".join(current))
    return sequences


def source_for(revision):
    if revision == "working-tree":
        return (ROOT / "src/main.cpp").read_text()
    return run_text(["git", "show", f"{revision}:src/main.cpp"])


def build_versions(build_dir, compiler, versions):
    build_dir.mkdir(parents=True, exist_ok=True)
    source_root = build_dir / "sources"
    binary_root = build_dir / "bin"
    source_root.mkdir(exist_ok=True)
    binary_root.mkdir(exist_ok=True)
    compiler_path = shutil.which(compiler)
    if compiler_path is None:
        raise RuntimeError(f"C++ compiler not found: {compiler}")

    metadata = []
    for label, revision, mode in versions:
        version_dir = source_root / label
        version_dir.mkdir(exist_ok=True)
        source_path = version_dir / "main.cpp"
        source_path.write_text(source_for(revision))
        binary_path = binary_root / f"sw_mapper-{label}"
        command = [compiler_path, *BUILD_FLAGS, str(source_path), "-o", str(binary_path)]
        print(f"BUILD {label} ({revision}, {mode})", flush=True)
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        (version_dir / "build.stdout.txt").write_text(completed.stdout)
        (version_dir / "build.stderr.txt").write_text(completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Build failed for {label}; see {version_dir / 'build.stderr.txt'}"
            )
        metadata.append((label, revision, mode, binary_path))

    (build_dir / "build-info.txt").write_text(
        "\n".join(
            [
                f"generated={datetime.datetime.now().astimezone().isoformat(timespec='seconds')}",
                f"platform={platform.platform()}",
                f"compiler={run_text([compiler_path, '--version']).splitlines()[0]}",
                f"flags={' '.join(BUILD_FLAGS)}",
                "historical versions use Git sources; v1.7.0-batch and final-auto use the working tree",
                "",
            ]
        )
    )
    return metadata


def parse_metrics(stderr, wall_seconds):
    app = re.search(
        r"\bmode=(\w+)(?: requested_mode=auto)? reads=(\d+) "
        r"reference_length=(\d+) seconds=([0-9.]+) "
        r"dp_seconds=([0-9.]+) GCUPS=([0-9.]+)",
        stderr,
    )
    time_line = re.search(
        r"(?m)^\s*([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys\s*$",
        stderr,
    )

    def integer(pattern):
        match = re.search(pattern, stderr, re.MULTILINE)
        return int(match.group(1)) if match else None

    if app is None:
        return {
            "mode": "",
            "reads": "",
            "reference_length": "",
            "program_seconds": None,
            "dp_seconds": None,
            "gcups": None,
            "user_seconds": None,
            "system_seconds": None,
            "max_rss_bytes": integer(r"([0-9]+)\s+maximum resident set size"),
            "page_faults": integer(r"([0-9]+)\s+page faults"),
            "page_reclaims": integer(r"([0-9]+)\s+page reclaims"),
            "instructions": integer(r"([0-9]+)\s+instructions retired"),
            "cycles": integer(r"([0-9]+)\s+cycles elapsed"),
            "ipc": None,
            "wall_seconds": wall_seconds,
            "cpu_cores": None,
        }

    user_seconds = float(time_line.group(2)) if time_line else None
    system_seconds = float(time_line.group(3)) if time_line else None
    cpu_seconds = (
        user_seconds + system_seconds
        if user_seconds is not None and system_seconds is not None else None
    )
    instructions = integer(r"([0-9]+)\s+instructions retired")
    cycles = integer(r"([0-9]+)\s+cycles elapsed")
    program_seconds = float(app.group(4))
    return {
        "mode": app.group(1),
        "reads": int(app.group(2)),
        "reference_length": int(app.group(3)),
        "program_seconds": program_seconds,
        "dp_seconds": float(app.group(5)),
        "gcups": float(app.group(6)),
        "user_seconds": user_seconds,
        "system_seconds": system_seconds,
        "max_rss_bytes": integer(r"([0-9]+)\s+maximum resident set size"),
        "page_faults": integer(r"([0-9]+)\s+page faults"),
        "page_reclaims": integer(r"([0-9]+)\s+page reclaims"),
        "instructions": instructions,
        "cycles": cycles,
        "ipc": instructions / cycles if instructions and cycles else None,
        "wall_seconds": wall_seconds,
        "cpu_cores": cpu_seconds / program_seconds
        if cpu_seconds is not None and program_seconds > 0 else None,
    }


def run_versions(metadata, reference, reads, output_dir, scoring, threads, repeats):
    result_root = output_dir / "results"
    result_root.mkdir(parents=True, exist_ok=True)
    (output_dir / "run-info.txt").write_text(
        "\n".join(
            [
                f"reference={reference}",
                f"reads={reads}",
                f"match={scoring[0]} mismatch={scoring[1]} gap={scoring[2]}",
                f"threads={threads} repeats={repeats}",
                "",
            ]
        )
    )
    rows = []
    baseline_output = None
    baseline_label = None

    for label, revision, mode, binary in metadata:
        version_result_dir = result_root / label
        version_result_dir.mkdir(exist_ok=True)
        for repeat in range(1, repeats + 1):
            command = [
                str(binary),
                "--reference", str(reference),
                "--reads", str(reads),
                "--mode", mode,
                "--match", str(scoring[0]),
                "--mismatch", str(scoring[1]),
                "--gap", str(scoring[2]),
            ]
            if mode in ("wavefront", "batch", "auto"):
                command.extend(("--threads", str(threads)))

            print(f"RUN {label} {repeat}/{repeats}: {shlex.join(command)}", flush=True)
            start = time.perf_counter()
            completed = subprocess.run(
                ["/usr/bin/time", "-l", *command],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            wall_seconds = time.perf_counter() - start
            run_stem = f"run-{repeat}"
            stdout_path = version_result_dir / f"{run_stem}.csv"
            stderr_path = version_result_dir / f"{run_stem}.time.txt"
            stdout_path.write_text(completed.stdout)
            stderr_path.write_text(completed.stderr)
            metric = parse_metrics(completed.stderr, wall_seconds)
            output_match = None
            if completed.returncode == 0:
                if baseline_output is None:
                    baseline_output = completed.stdout
                    baseline_label = label
                    output_match = True
                else:
                    output_match = completed.stdout == baseline_output
            row = {
                "version": label,
                "source": revision,
                "mode": mode,
                "repeat": repeat,
                "exit_code": completed.returncode,
                "output_matches_" + str(baseline_label): output_match,
                "stdout_csv": str(stdout_path.relative_to(result_root)),
                "stderr_time_log": str(stderr_path.relative_to(result_root)),
                **metric,
            }
            rows.append(row)
            if completed.returncode != 0:
                print(f"  ERROR: see {stderr_path}", flush=True)

    summary_path = result_root / "summary.csv"
    fields = (
        "version", "source", "mode", "repeat", "exit_code", "output_matches_baseline",
        "wall_seconds", "program_seconds", "dp_seconds", "gcups", "user_seconds",
        "system_seconds", "cpu_cores", "max_rss_bytes", "page_faults", "page_reclaims",
        "instructions", "cycles", "ipc", "stdout_csv", "stderr_time_log",
    )
    with summary_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            row["output_matches_baseline"] = row.pop("output_matches_" + str(baseline_label), None)
            writer.writerow(row)

    print("\nRESULT SUMMARY", flush=True)
    print("version                 mode       match  wall(s)   DP(s)   GCUPS   CPU(s)/cores   RSS(MiB)  IPC  CSV", flush=True)
    for label, _revision, mode, _binary in metadata:
        selected = [row for row in rows if row["version"] == label and row["exit_code"] == 0]
        if not selected:
            print(f"{label:23} {mode:10} BUILD/RUN ERROR", flush=True)
            continue
        numeric = ("wall_seconds", "program_seconds", "dp_seconds", "gcups", "user_seconds", "system_seconds", "cpu_cores", "max_rss_bytes", "ipc")
        medians = {}
        for key in numeric:
            values = [row[key] for row in selected if isinstance(row.get(key), (int, float))]
            medians[key] = statistics.median(values) if values else None
        cpu = (
            medians["user_seconds"] + medians["system_seconds"]
            if medians["user_seconds"] is not None and medians["system_seconds"] is not None else None
        )
        rss = medians["max_rss_bytes"] / (1024 * 1024) if medians["max_rss_bytes"] else None
        matches = all(row.get("output_matches_baseline") for row in selected)
        print(
            f"{label:23} {mode:10} "
            f"{'yes' if matches else 'NO':5} "
            f"{medians['wall_seconds'] or 0:7.3f} "
            f"{medians['dp_seconds'] or 0:7.3f} "
            f"{medians['gcups'] or 0:7.3f} "
            f"{(cpu or 0):6.2f}/{(medians['cpu_cores'] or 0):4.1f} "
            f"{(rss or 0):8.1f} "
            f"{(medians['ipc'] or 0):4.2f} "
            f"{selected[0]['stdout_csv']}",
            flush=True,
        )

    print(f"\nAll per-version CSV outputs and raw time logs: {result_root}", flush=True)
    print(f"Summary table: {summary_path}", flush=True)
    print(f"CSV comparison baseline: {baseline_label}", flush=True)
    print("First alignment result by version:", flush=True)
    for label, _revision, mode, _binary in metadata:
        selected = [row for row in rows if row["version"] == label and row["exit_code"] == 0]
        if not selected:
            continue
        first_csv = result_root / selected[0]["stdout_csv"]
        with first_csv.open(newline="") as input_file:
            first_row = next(csv.DictReader(input_file), None)
        if first_row:
            print(
                f"  {label:23} mode={mode:9} read={first_row['read_id']} "
                f"score={first_row['score']} ref=[{first_row['reference_start']},{first_row['reference_end']}) "
                f"read=[{first_row['read_start']},{first_row['read_end']})",
                flush=True,
            )
    mismatches = [row for row in rows if row["exit_code"] == 0 and row.get("output_matches_baseline") is False]
    if mismatches:
        print(f"WARNING: {len(mismatches)} run(s) differed from baseline {baseline_label}.", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build every Smith-Waterman iteration and run it on the same FASTA files."
    )
    parser.add_argument("--reference", type=pathlib.Path)
    parser.add_argument("--reads", type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, default=ROOT / "build/presentation-versions")
    parser.add_argument("--compiler", default="clang++")
    parser.add_argument("--threads", type=int, default=0, help="thread count for parallel versions (0=automatic)")
    parser.add_argument("--match", type=int, default=2)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-2)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--run-only", action="store_true", help="use binaries already built in --output-dir")
    parser.add_argument(
        "--versions",
        nargs="+",
        choices=tuple(label for label, _revision, _mode in VERSIONS),
        help="subset of version labels to build/run; default is all",
    )
    args = parser.parse_args()

    if platform.system() != "Darwin":
        parser.error("/usr/bin/time -l memory measurements require macOS")
    if args.build_only and args.run_only:
        parser.error("--build-only and --run-only cannot be used together")
    if args.repeats <= 0:
        parser.error("--repeats must be positive")
    if args.threads < 0:
        parser.error("--threads cannot be negative")
    if args.match <= 0 or args.mismatch > 0 or args.gap >= 0:
        parser.error("expected match > 0, mismatch <= 0, and gap < 0")

    selected_labels = set(args.versions or (label for label, _revision, _mode in VERSIONS))
    selected = [item for item in VERSIONS if item[0] in selected_labels]
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.build_only:
        if args.reference is None or args.reads is None:
            parser.error("--reference and --reads are required unless using --build-only")
        reference = args.reference.resolve()
        reads = args.reads.resolve()
        if not reference.is_file() or not reads.is_file():
            parser.error("reference and reads must point to existing FASTA files")
        if len(fasta_sequences(reference)) != 1:
            parser.error("reference FASTA must contain exactly one sequence")
        if not fasta_sequences(reads):
            parser.error("reads FASTA must contain at least one sequence")

    if args.run_only:
        metadata = []
        for label, revision, mode in selected:
            binary = output_dir / "bin" / f"sw_mapper-{label}"
            if not binary.is_file():
                parser.error(f"binary missing; build first: {binary}")
            metadata.append((label, revision, mode, binary))
    else:
        metadata = build_versions(output_dir, args.compiler, selected)

    if args.build_only:
        print(f"Built {len(metadata)} version binaries in {output_dir / 'bin'}")
        return 0

    run_versions(
        metadata,
        reference,
        reads,
        output_dir,
        (args.match, args.mismatch, args.gap),
        args.threads,
        args.repeats,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
