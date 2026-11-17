#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/run_version.sh VERSION REFERENCE.fa READS.fa

VERSION is one of:
  v1.0.0 v1.1.0 v1.2.0 v1.3.0 v1.4.0 v1.5.0 v1.6.0
  v1.7.0-batch final-auto

The script compiles that C++ source with Apple clang, runs it with the default
scoring (match=2, mismatch=-1, gap=-2), and records CSV plus /usr/bin/time -l
output under build/presentation-versions/results/VERSION/.

Set THREADS=N to override the automatic thread count in parallel modes.
EOF
}

if [[ $# -ne 3 ]]; then
  usage >&2
  exit 2
fi

version="$1"
reference="$2"
reads="$3"
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_root="$project_root/build/presentation-versions"
source_dir="$build_root/sources/$version"
binary_dir="$build_root/bin"
result_dir="$build_root/results/$version"
source_file="$source_dir/main.cpp"
binary="$binary_dir/sw_mapper-$version"
threads="${THREADS:-0}"

case "$version" in
  v1.0.0) revision="bc38b86"; mode="full" ;;
  v1.1.0) revision="bc38b86"; mode="rolling" ;;
  v1.2.0|v1.3.0|v1.4.0) revision="$version"; mode="trace" ;;
  v1.5.0|v1.6.0) revision="$version"; mode="wavefront" ;;
  v1.7.0-batch) revision="working-tree"; mode="batch" ;;
  final-auto) revision="working-tree"; mode="auto" ;;
  *) usage >&2; exit 2 ;;
esac

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This runner uses macOS /usr/bin/time -l." >&2
  exit 2
fi
if ! command -v clang++ >/dev/null 2>&1; then
  echo "clang++ is required to build the selected C++ version." >&2
  exit 2
fi
if [[ ! -f "$reference" || ! -f "$reads" ]]; then
  echo "Reference and reads must be existing FASTA files." >&2
  exit 2
fi

mkdir -p "$source_dir" "$binary_dir" "$result_dir"
if [[ "$revision" == "working-tree" ]]; then
  cp "$project_root/src/main.cpp" "$source_file"
else
  git -C "$project_root" show "$revision:src/main.cpp" > "$source_file"
fi

echo "Building $version from $revision ($mode)"
clang++ -O3 -DNDEBUG -g -std=c++17 -pthread "$source_file" -o "$binary"

command=(
  "$binary"
  --reference "$reference"
  --reads "$reads"
  --mode "$mode"
  --match 2
  --mismatch -1
  --gap -2
)
if [[ "$mode" == "wavefront" || "$mode" == "batch" || "$mode" == "auto" ]]; then
  command+=(--threads "$threads")
fi

echo "Running $version on the supplied FASTA files"
if /usr/bin/time -l "${command[@]}" \
  > "$result_dir/result.csv" 2> "$result_dir/time.txt"; then
  result_code=0
else
  result_code=$?
fi

echo "Alignment CSV: $result_dir/result.csv"
cat "$result_dir/result.csv"
echo "Timing and memory log: $result_dir/time.txt"
cat "$result_dir/time.txt"
exit "$result_code"
