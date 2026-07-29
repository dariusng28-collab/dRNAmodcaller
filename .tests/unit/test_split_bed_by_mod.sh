#!/usr/bin/env bash
set -euo pipefail

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cat > "$tmpdir/input.bed" <<'BED'
chr1	1	2	m	10	+	1	2	255,0,0	10	90	9	1	0	0	0	0	0
chr1	2	3	a	10	+	2	3	255,0,0	10	80	8	2	0	0	0	0	0
chr1	3	4	foo_m	10	+	3	4	255,0,0	10	70	7	3	0	0	0	0	0
BED

mods='[{"name":"m5C","code":"m","output_suffix":"m5C"},{"name":"m6A","code":"a","output_suffix":"m6A"}]'

python3 workflow/scripts/split_bed_by_mod.py \
  --input-bed "$tmpdir/input.bed" \
  --modifications-json "$mods" \
  --output-dir "$tmpdir/out" \
  --sample sample1 \
  --log "$tmpdir/split.log"

test "$(wc -l < "$tmpdir/out/sample1.m5C.filtered.bed")" -eq 1
test "$(wc -l < "$tmpdir/out/sample1.m6A.filtered.bed")" -eq 1
! grep -q 'foo_m' "$tmpdir/out/sample1.m5C.filtered.bed"
