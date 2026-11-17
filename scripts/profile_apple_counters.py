#!/usr/bin/env python3
"""Capture macOS powermetrics alongside /usr/bin/time -l for one command.

Run with sudo so powermetrics can sample per-process CPU statistics. The target
program is dropped back to the invoking user when this script is launched via
sudo, so benchmark input/output files keep normal user permissions.
"""

import argparse
import datetime
import os
import pathlib
import platform
import pwd
import signal
import subprocess
import sys
import time


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Collect Apple powermetrics task/CPU samples and /usr/bin/time -l "
            "for a command. Invoke with sudo; the target runs as the caller."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        help="output directory (default: build/performance-counters/<timestamp>)",
    )
    parser.add_argument(
        "--sample-rate-ms",
        type=int,
        default=50,
        help="powermetrics sample interval in milliseconds (default: 50)",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("provide the command after --")
    if args.sample_rate_ms <= 0:
        parser.error("--sample-rate-ms must be positive")
    if os.geteuid() != 0:
        parser.error("run through sudo so powermetrics can access task counters")

    try:
        caller_uid = int(os.environ["SUDO_UID"])
        caller_gid = int(os.environ["SUDO_GID"])
    except (KeyError, ValueError):
        parser.error("run as `sudo python3 scripts/profile_apple_counters.py ...`")
    caller = pwd.getpwuid(caller_uid)

    output_dir = args.output_dir
    if output_dir is None:
        stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        output_dir = pathlib.Path("build/performance-counters") / stamp
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    power_log = output_dir / "powermetrics.txt"
    power_error_log = output_dir / "powermetrics-errors.txt"
    time_log = output_dir / "time.txt"
    stdout_log = output_dir / "stdout.txt"
    metadata_log = output_dir / "metadata.txt"
    sample_command = [
        "/usr/bin/powermetrics",
        "--samplers", "tasks,cpu_power",
        "--show-process-ipc",
        "--sample-rate", str(args.sample_rate_ms),
        "--sample-count", "-1",
        "--buffer-size", "1",
        "--output-file", str(power_log),
    ]

    def return_to_caller() -> None:
        os.initgroups(caller.pw_name, caller_gid)
        os.setgid(caller_gid)
        os.setuid(caller_uid)

    power_sampler = None
    target = None
    started = time.monotonic()
    return_code = 1
    try:
        with power_error_log.open("w") as power_errors:
            power_sampler = subprocess.Popen(
                sample_command,
                stdout=subprocess.DEVNULL,
                stderr=power_errors,
            )
        with stdout_log.open("w") as stdout_file, time_log.open("w") as time_file:
            target = subprocess.Popen(
                ["/usr/bin/time", "-l", *command],
                stdout=stdout_file,
                stderr=time_file,
                preexec_fn=return_to_caller,
            )
            return_code = target.wait()
    except KeyboardInterrupt:
        if target is not None and target.poll() is None:
            target.send_signal(signal.SIGINT)
            target.wait()
        return_code = 130
    finally:
        if power_sampler is not None and power_sampler.poll() is None:
            power_sampler.send_signal(signal.SIGINT)
            try:
                power_sampler.wait(timeout=10)
            except subprocess.TimeoutExpired:
                power_sampler.kill()
                power_sampler.wait()

    elapsed = time.monotonic() - started
    metadata_log.write_text(
        "\n".join(
            [
                f"command={' '.join(command)}",
                f"caller={caller.pw_name} uid={caller_uid} gid={caller_gid}",
                f"platform={platform.platform()}",
                f"architecture={platform.machine()}",
                f"sample_rate_ms={args.sample_rate_ms}",
                f"elapsed_seconds={elapsed:.6f}",
                f"exit_code={return_code}",
                "powermetrics=tasks,cpu_power,per-process IPC",
                "time=/usr/bin/time -l (includes instructions and cycles on Apple Silicon)",
                "",
            ]
        )
    )

    for path in output_dir.iterdir():
        os.chown(path, caller_uid, caller_gid)
    os.chown(output_dir, caller_uid, caller_gid)
    print(f"output={output_dir}")
    print(f"exit_code={return_code}")
    print("powermetrics includes per-process CPU/IPC and frequency samples.")
    print("It does not expose cache-miss or memory-stall counters.")
    return return_code


if __name__ == "__main__":
    sys.exit(main())
