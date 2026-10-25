#!/usr/bin/env python3
"""Keep the machine at a number of SCycle jobs, starting the runs of several batches in order.

Every --interval seconds this counts the SCycle processes on the machine (an executable named main
given an input *.in, inside or outside the batches; each MPI rank counts once). While fewer than
--jobs run, it gives one free slot to the first batch, in the order given, with a run waiting for
one (new, stopped or queued, as tools/batch.sh status reports): it starts that batch's queue with
one job, or, when the queue runs already, restarts it with one job more than the batch runs
(tools/batch.sh stop, then start: the running runs resume from their checkpoints). A slot counts as
free only if it is still free --settle seconds later, since a run just started, or a queue moving on
to its next run, takes a moment to show. Failed runs are left alone. Every batch must have its
executable already (<dir>/bin/main, copied by a first tools/batch.sh start with SCYCLE_BIN), so that
a batch never takes whatever source/main is. It exits when no run is left waiting. Run it detached.

Usage: fill_slots.py BATCH [BATCH ...] [--jobs 4] [--interval 300] [--settle 60] [--dry-run]
"""
import argparse
import datetime
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
BATCH_SH = os.path.join(ROOT, 'tools', 'batch.sh')


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def scycle_processes():
    """SCycle processes running on the machine."""
    out = subprocess.run(['ps', '-axo', 'command='], capture_output=True, text=True).stdout
    return sum(1 for w in (line.split() for line in out.splitlines())
               if len(w) >= 2 and os.path.basename(w[0]) == 'main' and w[1].endswith('.in'))


def status(batch):
    """Whether the batch's queue runs, and the state of each run, from tools/batch.sh status."""
    lines = subprocess.run([BATCH_SH, 'status', batch], capture_output=True, text=True, cwd=ROOT).stdout.splitlines()
    states = {w[0]: w[1] for w in (line.split() for line in lines[2:]) if len(w) >= 2}
    return bool(lines) and lines[0] == 'queue: running', states


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('batches', nargs='+')
    a.add_argument('--jobs', type=int, default=4, help='SCycle processes to keep running')
    a.add_argument('--interval', type=float, default=300.0, help='seconds between checks')
    a.add_argument('--settle', type=float, default=60.0, help='seconds a slot must stay free')
    a.add_argument('--dry-run', action='store_true', help='report once, change nothing')
    o = a.parse_args()
    batches = [os.path.abspath(b) for b in o.batches]
    missing = [b for b in batches if not os.access(os.path.join(b, 'bin', 'main'), os.X_OK)]
    if missing:
        sys.exit('no executable in ' + ', '.join(os.path.join(os.path.relpath(b, ROOT), 'bin', 'main') for b in missing)
                 + ': start each batch once with SCYCLE_BIN=<build> tools/batch.sh start, or copy one there')
    print(f'{now()} keeping {o.jobs} SCycle jobs; batches in order: '
          + ', '.join(os.path.relpath(b, ROOT) for b in batches), flush=True)
    while True:
        n, waiting = scycle_processes(), []
        if n < o.jobs and not o.dry_run:
            time.sleep(o.settle)
            n = scycle_processes()
        for b in batches:
            queue, st = status(b)
            w = [r for r, s in st.items() if s in ('new', 'stopped', 'waiting')]
            if w:
                waiting.append((b, queue, sum(s == 'running' for s in st.values()), w))
        if o.dry_run:
            print(f'  {n} SCycle processes; waiting: '
                  + ('; '.join(f'{os.path.relpath(b, ROOT)}: {" ".join(w)}' for b, _, _, w in waiting) or 'none'))
            return
        if not waiting:
            print(f'{now()} no run is left waiting; done', flush=True)
            return
        if n < o.jobs:
            b, queue, running, w = waiting[0]
            jobs = running + 1 if queue else 1
            print(f'{now()} {n} SCycle jobs: {"restarting" if queue else "starting"} {os.path.relpath(b, ROOT)} '
                  f'with {jobs} job(s) (waiting: {" ".join(w)})', flush=True)
            if queue:
                subprocess.run([BATCH_SH, 'stop', b], cwd=ROOT)
            subprocess.run([BATCH_SH, 'start', b, '-j', str(jobs)], cwd=ROOT)
            sys.stdout.flush()
        time.sleep(o.interval)


if __name__ == '__main__':
    sys.exit(main())
