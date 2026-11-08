#!/usr/bin/env python3
"""Generate the deterministic DNA dataset used by the README benchmark."""

import argparse
import pathlib
import random


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    parser.add_argument("--reference-length", type=int, default=12_000)
    parser.add_argument("--read-count", type=int, default=64)
    parser.add_argument("--read-length", type=int, default=160)
    parser.add_argument("--min-read-length", type=int)
    parser.add_argument("--max-read-length", type=int)
    parser.add_argument("--substitutions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=21021)
    args = parser.parse_args()

    if (args.min_read_length is None) != (args.max_read_length is None):
        parser.error("provide both min-read-length and max-read-length, or neither")
    min_length = (
        args.read_length if args.min_read_length is None else args.min_read_length
    )
    max_length = (
        args.read_length if args.max_read_length is None else args.max_read_length
    )
    if args.reference_length <= max_length:
        parser.error("reference length must be greater than the longest read")
    if args.read_count <= 0 or min_length <= 0 or max_length < min_length:
        parser.error("read count and read lengths must be positive and ordered")
    if args.substitutions < 0 or args.substitutions > min_length:
        parser.error("substitutions must be between zero and the shortest read")

    rng = random.Random(args.seed)
    alphabet = "ACGT"
    reference = "".join(rng.choice(alphabet) for _ in range(args.reference_length))
    reads = []
    # The exclusive upper bound intentionally matches the first recorded run.
    for index in range(args.read_count):
        read_length = (
            min_length
            if min_length == max_length
            else rng.randint(min_length, max_length)
        )
        start = rng.randrange(0, len(reference) - read_length)
        bases = list(reference[start : start + read_length])
        for position in rng.sample(range(len(bases)), args.substitutions):
            bases[position] = rng.choice([base for base in alphabet if base != bases[position]])
        reads.append((f"read_{index:03d}", "".join(bases)))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "reference.fa").write_text(f">reference\n{reference}\n")
    (args.output_dir / "reads.fa").write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in reads)
    )
    print(
        f"reference={args.reference_length}bp reads={args.read_count} "
        f"read_lengths={min_length}-{max_length}bp "
        f"substitutions={args.substitutions} seed={args.seed}"
    )


if __name__ == "__main__":
    main()
