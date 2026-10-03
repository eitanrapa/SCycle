#!/usr/bin/env bash
# Run a batch of SCycle inputs detached from the terminal, so the runs outlive the shell (or ssh
# session) that started them, and resume them after a stop or a reboot.
#
# A batch is a directory of input files <dir>/<run>.in (examples/two_faults/stage5_batch1.py writes
# one). Each run writes to the outputDir its input names, and its log to <outputDir>run.log; relative
# outputDir paths are taken from the repository root.
#
# Usage:
#   tools/batch.sh start  <dir> [-j N] [-n R] [--wait]
#       start every run that is neither running nor done, at most N at a time (default: all at once),
#       each on R MPI ranks (default 1). Run it again after a stop or a reboot: each run continues
#       from its last checkpoint. --wait stays in the foreground until the runs end (for the job
#       script of a scheduler such as Slurm).
#   tools/batch.sh status <dir>     state, step and simulated time of each run
#   tools/batch.sh stop   <dir>     stop the queue and the runs; start resumes them
#
# stop sends SIGTERM to SCycle, which finishes its step, writes a checkpoint, closes its output and
# exits with status 143; a batch system's SIGTERM (Slurm at the time limit, a shutdown) does the same.
# Do not kill -9 a run: HDF5 output killed while being written can become unreadable.
# The first start copies the executable (SCYCLE_BIN, default source/main; build it with an optimized
# PETSC_ARCH) to <dir>/bin/main and records where it came from in <dir>/bin/BUILD.txt. Later starts
# use that copy, so rebuilding the code does not change a batch in progress.
# A run that ends on its own leaves its exit status in <outputDir>run.exit: 0 is finished, anything
# else failed (the reason is in run.log). start skips both; delete run.exit to run one again (raise
# maxTime in its input first to extend a finished run). Resuming needs restartFromChkpt = 1, the
# default.
#
# Environment:
#   SCYCLE_BIN     executable copied by the first start (default: source/main)
#   MPIEXEC        MPI launcher for R > 1 (default: mpirun)
#   MPIEXEC_FLAGS  extra launcher flags; with Open MPI and several multi-rank runs at once, use
#                  "--bind-to none", or each run's ranks are bound to the same first cores
#   OMP_NUM_THREADS, OPENBLAS_NUM_THREADS   default 1, so a threaded BLAS does not oversubscribe

self=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
root=$(cd "$(dirname "$0")/.." && pwd)

usage() { sed -n '/^# Usage:/,/^#$/p' "$self" | sed 's/^# \{0,1\}//' >&2; exit 2; }

name_of() { basename "$1" .in; }

# the outputDir (a file-name prefix) that input $1 names
prefix_of() {
  local p
  p=$(sed -n 's/^outputDir = \([^ ]*\).*/\1/p' "$1" | tail -1)
  [[ -n $p ]] || { echo "no outputDir in $1" >&2; return 1; }
  [[ $p == /* ]] || p=$root/$p
  printf '%s\n' "$p"
}

# true if file $1 holds the pid of a live process whose command line contains $2
# (a pid left by a reboot may belong to another process by now)
alive() {
  local pid cmd
  pid=$(cat "$1" 2> /dev/null) && [[ -n $pid ]] || return 1
  cmd=$(ps -p "$pid" -o command= 2> /dev/null) && [[ $cmd == *"$2"* ]]
}

# true while this batch's queue runs: the pid in queue.pid is first nohup, then setsid (or perl), then
# bash running _queue, then xargs, whose command lines all hold "batch.sh _queue|_run <dir> "
queue_alive() {
  local pid cmd
  pid=$(cat "$dir/queue.pid" 2> /dev/null) && [[ -n $pid ]] || return 1
  cmd=$(ps -p "$pid" -o command= 2> /dev/null) && [[ $cmd == *"batch.sh _"*" $dir "* ]]
}

# $1 with the characters special in an extended regular expression escaped
regex() { sed 's/[][\.*^$+?(){}|]/\\&/g' <<< "$1"; }

# running, finished, failed, waiting (queued), stopped (started before) or new
state_of() {
  local in=$1 p
  p=$(prefix_of "$in") || { echo invalid; return; }
  if alive "${p}run.pid" "$in"; then echo running
  elif [[ -f ${p}run.exit ]]; then [[ $(cat "${p}run.exit") == 0 ]] && echo finished || echo failed
  elif queue_alive && grep -q -x -F -- "$in" "$dir/queue.txt" 2> /dev/null; then echo waiting
  elif [[ -f ${p}run.log ]]; then echo stopped
  else echo new; fi
}

start() {
  local jobs=0 np=1 wait=0 in p s pending=()
  while [[ $# -gt 0 ]]; do
    case $1 in
      -j) jobs=$2; shift 2 ;;
      -n) np=$2; shift 2 ;;
      --wait) wait=1; shift ;;
      *) usage ;;
    esac
  done
  [[ $jobs =~ ^[0-9]+$ && $np =~ ^[1-9][0-9]*$ ]] || usage
  if queue_alive; then echo "a queue is already running for $dir; stop it first" >&2; exit 1; fi

  mkdir -p "$dir/bin"
  if [[ ! -x $dir/bin/main ]]; then
    local src=${SCYCLE_BIN:-$root/source/main}
    [[ -x $src ]] || { echo "missing executable $src (run: make -C source)" >&2; exit 2; }
    cp "$src" "$dir/bin/main" || exit 2
    { echo "copied from $src on $(date '+%Y-%m-%d %H:%M') ($(uname -n))"
      echo "source: $(git -C "$root" rev-parse --abbrev-ref HEAD 2> /dev/null)" \
           "$(git -C "$root" describe --always --dirty 2> /dev/null)"
      { ldd "$src" || otool -L "$src"; } 2> /dev/null | grep -m1 -o '/[^ ]*libpetsc[^ ]*'
    } > "$dir/bin/BUILD.txt"
  fi
  echo "executable: $dir/bin/main, $(head -1 "$dir/bin/BUILD.txt" 2> /dev/null)"

  for in in "${inputs[@]}"; do
    s=$(state_of "$in")
    case $s in
      new) pending+=("$in"); echo "$(name_of "$in"): start" ;;
      stopped) pending+=("$in"); echo "$(name_of "$in"): resume" ;;
      failed) p=$(prefix_of "$in")
              echo "$(name_of "$in"): failed (exit $(cat "${p}run.exit")), see ${p}run.log; delete ${p}run.exit to retry" ;;
      *) echo "$(name_of "$in"): $s" ;;
    esac
  done
  [[ ${#pending[@]} -gt 0 ]] || { echo "nothing to start"; return 0; }
  [[ $jobs -gt 0 ]] || jobs=${#pending[@]}
  local cores; cores=$(getconf _NPROCESSORS_ONLN 2> /dev/null || echo 0)
  (( cores > 0 && jobs * np > cores )) && echo "note: up to $((jobs * np)) processes on $cores cores"

  printf '%s\n' "${pending[@]}" > "$dir/queue.txt"
  if (( wait )); then
    bash "$self" _queue "$dir" "$jobs" "$np" < /dev/null
    status
  else
    # a new session that ignores hangups outlives this shell; each exec keeps the pid, which is
    # recorded here so that a start right after this one already sees the queue
    if command -v setsid > /dev/null; then
      nohup setsid bash "$self" _queue "$dir" "$jobs" "$np" < /dev/null >> "$dir/queue.log" 2>&1 &
    else   # macOS has no setsid
      nohup perl -MPOSIX -e 'POSIX::setsid(); exec @ARGV' bash "$self" _queue "$dir" "$jobs" "$np" \
        < /dev/null >> "$dir/queue.log" 2>&1 &
    fi
    echo $! > "$dir/queue.pid"
    echo "queue: ${#pending[@]} runs, $jobs at a time, $np rank(s) each; follow it with: $0 status $dir"
  fi
}

status() {
  local in p s last step t tmax idle
  if queue_alive; then echo "queue: running"; else echo "queue: none"; fi
  printf '%-12s %-10s %10s %10s %10s\n' run state step "t (yr)" "of (yr)"
  for in in "${inputs[@]}"; do
    p=$(prefix_of "$in") || continue
    s=$(state_of "$in")
    [[ $s == failed ]] && s="failed($(cat "${p}run.exit"))"
    idle=""
    [[ $s == running && -n $(find "${p}run.log" -mmin +30 2> /dev/null) ]] && idle="  no output for 30 min"
    last=$(tail -c 20000 "${p}run.log" 2> /dev/null | grep -E '^[0-9]+: t = ' | tail -1)
    step=${last%%:*}
    t=$(sed -n 's/^[0-9]*: t = \([^ ,]*\) s.*/\1/p' <<< "$last")
    tmax=$(sed -n 's/^maxTime = \([^ ]*\).*/\1/p' "$in" | tail -1)
    awk -v n="$(name_of "$in")" -v s="$s" -v k="${step:--}" -v t="$t" -v m="$tmax" -v x="$idle" 'BEGIN {
      Y = 3.15576e7
      printf "%-12s %-10s %10s %10s %10s%s\n", n, s, k, (t == "" ? "-" : sprintf("%.1f", t/Y)),
             (m == "" ? "-" : sprintf("%.0f", m/Y)), x }'
  done
}

stop() {
  local in p msg left waited=0 running=()
  if queue_alive; then kill "$(cat "$dir/queue.pid")" && echo "queue: stopped, no new runs start"; fi
  for in in "${inputs[@]}"; do
    p=$(prefix_of "$in") || continue
    alive "${p}run.pid" "$in" || continue
    touch "${p}run.stopping"   # tells the run's wrapper that this exit is a stop, not a failure
    # signal SCycle itself, every rank of it, rather than an MPI launcher, which may kill the ranks
    # before they have closed their files
    pkill -TERM -f "^$(regex "$dir/bin/main") $(regex "$in")\$" || kill "$(cat "${p}run.pid")"
    running+=("$in")
  done
  [[ ${#running[@]} -gt 0 ]] || { echo "no run is running"; return 0; }
  echo "asked ${#running[@]} run(s) to stop; each finishes its step and writes a checkpoint"
  while :; do
    left=()
    for in in "${running[@]}"; do
      p=$(prefix_of "$in")
      alive "${p}run.pid" "$in" && left+=("$in")
    done
    [[ ${#left[@]} -gt 0 && $waited -lt 600 ]] || break
    sleep 2; waited=$((waited + 2))
  done
  for in in "${running[@]}"; do
    p=$(prefix_of "$in")
    if alive "${p}run.pid" "$in"; then echo "$(name_of "$in"): still running after $waited s"; continue; fi
    # SCycle's report of the stop, from the part of the log written since the run's last start
    msg=$(awk '/^== .* start on /{m=""} /^Stopping on signal/{m=$0} END{print m}' "${p}run.log")
    if [[ -n $msg ]]; then echo "$(name_of "$in"): $(sed 's/^Stopping on signal [0-9]* at /stopped at /' <<< "$msg")"
    else echo "$(name_of "$in"): stopped while starting, before its first step"; fi
  done
}

# internal: the queue, at most $1 runs at a time, each on $2 ranks
queue() {
  echo $$ > "$dir/queue.pid"
  exec xargs -P "$1" -I{} bash "$self" _run "$dir" "$2" {} < "$dir/queue.txt"
}

# internal: one run of input $2 on $1 ranks
run() {
  local np=$1 in=$2 p pid rc off
  p=$(prefix_of "$in") || return 0
  alive "${p}run.pid" "$in" && return 0      # started meanwhile by another queue
  [[ -f ${p}run.exit ]] && return 0
  mkdir -p "$(dirname "${p}x")"
  rm -f "${p}run.stopping"
  echo "== $(date '+%Y-%m-%d %H:%M:%S') start on $(uname -n): $dir/bin/main, $np rank(s)" >> "${p}run.log"
  off=$(wc -c < "${p}run.log")   # where this start's part of the log begins
  export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
  cd "$root" || return 0
  trap ':' TERM INT   # outlive a signal meant for the whole job, to record how the run ended
  if (( np > 1 )); then
    ${MPIEXEC:-mpirun} ${MPIEXEC_FLAGS:-} -n "$np" "$dir/bin/main" "$in" >> "${p}run.log" 2>&1 < /dev/null &
  else
    "$dir/bin/main" "$in" >> "${p}run.log" 2>&1 < /dev/null &
  fi
  pid=$!
  echo "$pid" > "${p}run.pid"
  while :; do
    wait "$pid"; rc=$?
    kill -0 "$pid" 2> /dev/null || break    # a trapped signal ends wait early while the run goes on
  done
  # stopped by this script, by a signal SCycle caught (it reports the stop), or killed by one
  if [[ -f ${p}run.stopping ]] || tail -c +$((off + 1)) "${p}run.log" | grep -q '^Stopping on signal' ||
     [[ $rc == 129 || $rc == 130 || $rc == 137 || $rc == 143 ]]; then
    rm -f "${p}run.stopping"
    echo "== $(date '+%Y-%m-%d %H:%M:%S') stopped (exit $rc)" >> "${p}run.log"
  else
    echo "$rc" > "${p}run.exit"
    echo "== $(date '+%Y-%m-%d %H:%M:%S') exit $rc" >> "${p}run.log"
  fi
}

cmd=${1:-}
[[ -n ${2:-} && -d ${2:-} ]] || usage
dir=$(cd "$2" && pwd)
shift 2
inputs=("$dir"/*.in)
[[ -f ${inputs[0]} ]] || { echo "no input files (*.in) in $dir" >&2; exit 2; }

case $cmd in
  start) start "$@" ;;
  status) status ;;
  stop) stop ;;
  _queue) queue "$@" ;;
  _run) run "$@" ;;
  *) usage ;;
esac
