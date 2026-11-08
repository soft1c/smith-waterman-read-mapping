#!/usr/bin/env python3
"""Record peak RSS for each implementation on every benchmark profile (macOS)."""

import argparse
import csv
import pathlib
import re
import subprocess


def run_peak_rss(command):
    process = subprocess.run(
        ["/usr/bin/time", "-l", *map(str, command)],
        capture_output=True,
        text=True,
        check=True,
    )
    match = re.search(r"([0-9]+)\s+maximum resident set size", process.stderr)
    if not match:
        raise RuntimeError(f"Could not parse macOS peak RSS output: {process.stderr}")
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--suite-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    args = parser.parse_args()

    rows = []
    for profile in ("small", "medium", "large"):
        profile_dir = args.suite_dir / profile
        data_dir = profile_dir / "data"
        online_binary = profile_dir / "online-build/online_sw_sequential"
        pair_file = profile_dir / "online-build/online-input.txt"
        for mode in ("online", "full", "rolling", "simd"):
            if mode == "online":
                command = [
                    online_binary,
                    "-path",
                    pair_file,
                    "-id",
                    "memory",
                    "-match",
                    "2",
                    "-mismatch",
                    "-1",
                    "-gap",
                    "-2",
                ]
            else:
                command = [
                    args.binary,
                    "--reference",
                    data_dir / "reference.fa",
                    "--reads",
                    data_dir / "reads.fa",
                    "--mode",
                    mode,
                ]
            peak_bytes = run_peak_rss(command)
            row = {
                "profile": profile,
                "mode": mode,
                "peak_rss_bytes": peak_bytes,
                "peak_rss_mib": peak_bytes / (1024 * 1024),
                "samples": 1,
            }
            rows.append(row)
            print(
                f"profile={profile} mode={mode} "
                f"peak_RSS={row['peak_rss_mib']:.2f} MiB"
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=("profile", "mode", "peak_rss_bytes", "peak_rss_mib", "samples"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"memory_csv={args.output.resolve()}")


if __name__ == "__main__":
    main()
