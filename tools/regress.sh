#!/usr/bin/env bash
# Regression check for SCycle.
# Runs example inputs into a scratch directory and compares their HDF5 output
# against a stored baseline with h5diff.
#
# Usage:
#   tools/regress.sh baseline <dir>   run the cases and store their output in <dir>
#   tools/regress.sh compare  <dir>   run the cases and compare against <dir>
#
# Environment:
#   SCYCLE_BIN     executable (default: source/main)
#   REGRESS_NP     MPI ranks (default 1; >1 runs through mpirun)
#   REGRESS_CASES  example names under examples/, or paths to other input files
#                  (default "ex1 ex2"); a path's case name is its file name without .in
#   REGRESS_DELTA  absolute tolerance for h5diff --delta (default: exact)
#   REGRESS_WORK   scratch directory for the runs (default: mktemp -d)
#
# Each case always starts fresh: outputDir is redirected to the scratch
# directory and checkpoint restarts are disabled, whatever the example says.
# checkpoint.h5 is not compared: it only holds the state of the last checkpoint.
set -uo pipefail

mode=${1:-}; ref=${2:-}
if [[ "$mode" != baseline && "$mode" != compare ]] || [[ -z "$ref" ]]; then
  echo "usage: $0 baseline|compare <dir>" >&2; exit 2
fi
root=$(cd "$(dirname "$0")/.." && pwd)
bin=${SCYCLE_BIN:-$root/source/main}
np=${REGRESS_NP:-1}
cases=${REGRESS_CASES:-"ex1 ex2"}
work=${REGRESS_WORK:-$(mktemp -d)}

[[ -x "$bin" ]] || { echo "missing executable $bin (run: make -C source)" >&2; exit 2; }
command -v h5diff >/dev/null || { echo "h5diff not found on PATH" >&2; exit 2; }
[[ "$mode" == compare && ! -d "$ref" ]] && { echo "no baseline at $ref" >&2; exit 2; }

status=0
for case in $cases; do
  if [[ -f "$case" ]]; then src=$case; ex=$(basename "$case" .in)
  else src="$root/examples/$case.in"; ex=$case; fi
  [[ -f "$src" ]] || { echo "FAIL $case: no input file $src"; status=1; continue; }
  out="$work/$ex"; rm -rf "$out"; mkdir -p "$out"
  in="$work/$ex.in"
  sed -e "s|^outputDir = .*|outputDir = $out/|" \
      -e '/^restartFromChkpt/d' -e '/^restartFromChkptSS/d' -e '/^retartFromChkptSS/d' \
      "$src" > "$in"
  printf 'restartFromChkpt = 0\nrestartFromChkptSS = 0\n' >> "$in"

  start=$SECONDS
  if (( np > 1 )); then mpirun -n "$np" "$bin" "$in" > "$out/stdout.log" 2>&1
  else "$bin" "$in" > "$out/stdout.log" 2>&1; fi
  rc=$?
  if (( rc != 0 )) || grep -q "PETSC ERROR\|Assertion failed" "$out/stdout.log"; then
    echo "FAIL $ex: run failed (exit $rc), see $out/stdout.log"; status=1; continue
  fi
  echo "ran  $ex in $((SECONDS - start)) s"

  if [[ "$mode" == baseline ]]; then
    mkdir -p "$ref/$ex"
    find "$out" -name '*.h5' ! -name checkpoint.h5 -exec cp {} "$ref/$ex/" \;
    cp "$out"/*.txt "$out/stdout.log" "$ref/$ex/" 2>/dev/null
    echo "     baseline stored in $ref/$ex"
    continue
  fi

  for f in "$ref/$ex"/*.h5; do
    name=$(basename "$f")
    if [[ ! -f "$out/$name" ]]; then echo "FAIL $ex/$name: missing from new run"; status=1; continue; fi
    # h5diff exits 0 even when datasets have different shapes (e.g. a run that stopped early);
    # it only prints "not comparable", so treat that as a difference too
    report=$(h5diff -r -c ${REGRESS_DELTA:+--delta=$REGRESS_DELTA} "$f" "$out/$name" 2>&1)
    rc=$?
    if (( rc == 0 )) && grep -qi "not comparable" <<< "$report"; then rc=1; fi
    case $rc in
      0) echo "ok   $ex/$name" ;;
      1) echo "DIFF $ex/$name"; grep -E "differences found|not comparable|dimensions" <<< "$report" | head -20; status=1 ;;
      *) echo "ERR  $ex/$name"; head -5 <<< "$report"; status=1 ;;
    esac
  done
done
echo "run directory: $work"
exit $status
