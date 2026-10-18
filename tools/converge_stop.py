#!/usr/bin/env python3
"""Stop the runs of a batch once their earthquake cycle has converged, and widen the batch later.

A spin-up from a steady state approaches its limit cycle slowly: the recurrence interval grows or
shrinks a little from one earthquake to the next. This watches the runs of a tools/batch.sh batch
and stops each one (SIGTERM: it writes a checkpoint) once its last --cycles recurrence intervals
each changed by less than --tol years; it then writes run.exit = 0, so the batch treats the run as
finished and its queue moves on, and notes the stop in the run's run.log. A run reaching its own
maxTime ends as usual.

With --widen-after OTHER --jobs N, once every run of the batch OTHER has finished it restarts this
batch with N jobs at a time (tools/batch.sh stop, then start: running runs resume from checkpoints).

Earthquakes are read from run.log: the first step shorter than 1 s after one longer than 1e5 s
(the step falls from interseismic sizes to coseismic ones). Run detached; it exits when no watched
run is left to stop and the widening, if any, is done.

Usage: converge_stop.py BATCH [--runs NAME ...] [--tol 0.5] [--cycles 2] [--interval 300]
                        [--widen-after OTHER --jobs N] [--dry-run]
"""
import argparse
import datetime
import glob
import os
import re
import signal
import subprocess
import sys
import time

YEAR = 3.15576e7
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
STEP = re.compile(r'^(\d+): t = (\S+) s, dt = (\S+),')


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def output_dir(inp):
    """The outputDir the input names (relative paths are taken from the repository root)."""
    p = None
    for line in open(inp):
        m = re.match(r'^outputDir = (\S+)', line)
        if m:
            p = m.group(1)
    if p is None:
        return None
    return p if os.path.isabs(p) else os.path.join(ROOT, p)


def earthquakes(log):
    """Onset times (s) of the earthquakes in a run.log."""
    onsets, big = [], False
    with open(log, errors='replace') as f:
        for line in f:
            m = STEP.match(line)
            if not m:
                continue
            t, dt = float(m.group(2)), float(m.group(3))
            if dt > 1e5:
                big = True
            elif dt < 1.0 and big:
                onsets.append(t)
                big = False
    return onsets


def alive(pidfile, inp):
    """The pid in pidfile if it is a live process whose command line names the input, else None."""
    try:
        pid = int(open(pidfile).read().strip())
        cmd = subprocess.run(['ps', '-p', str(pid), '-o', 'command='], capture_output=True, text=True).stdout
    except (OSError, ValueError):
        return None
    return pid if inp in cmd else None


def batch_finished(batch):
    """True when every input of the batch has run.exit 0 and none is running."""
    for inp in sorted(glob.glob(os.path.join(batch, '*.in'))):
        p = output_dir(inp)
        try:
            if open(p + 'run.exit').read().strip() != '0':
                return False
        except OSError:
            return False
        if alive(p + 'run.pid', inp):
            return False
    return True


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('batch')
    a.add_argument('--runs', nargs='*', help='run names to watch (default: every input of the batch)')
    a.add_argument('--tol', type=float, default=0.5, help='largest change of the recurrence interval (yr)')
    a.add_argument('--cycles', type=int, default=2, help='consecutive cycles within --tol')
    a.add_argument('--interval', type=float, default=300.0, help='seconds between checks')
    a.add_argument('--widen-after', help='batch whose end lets this one take more jobs')
    a.add_argument('--jobs', type=int, help='jobs at a time after --widen-after has finished')
    a.add_argument('--dry-run', action='store_true', help='report once, change nothing')
    o = a.parse_args()
    batch = os.path.abspath(o.batch)
    inputs = sorted(glob.glob(os.path.join(batch, '*.in')))
    if o.runs:
        inputs = [i for i in inputs if os.path.basename(i)[:-3] in o.runs]
    widened = o.widen_after is None
    print(f'{now()} watching {len(inputs)} runs of {batch}: tolerance {o.tol} yr over {o.cycles} cycles'
          + ('' if widened else f'; {o.jobs} jobs once {o.widen_after} has finished'), flush=True)
    while True:
        left = 0
        for inp in inputs:
            name, p = os.path.basename(inp)[:-3], output_dir(inp)
            if os.path.exists(p + 'run.exit'):
                continue                                  # finished, failed or stopped here
            left += 1
            if not os.path.exists(p + 'run.log'):
                continue                                  # not started yet
            ev = earthquakes(p + 'run.log')
            rec = [(ev[i + 1] - ev[i]) / YEAR for i in range(len(ev) - 1)]
            ch = [abs(rec[i + 1] - rec[i]) for i in range(len(rec) - 1)]
            done = len(ch) >= o.cycles and max(ch[-o.cycles:]) < o.tol
            if o.dry_run:
                print(f'  {name}: {len(ev)} earthquakes, recurrence {[round(r, 1) for r in rec[-4:]]} yr, '
                      f'changes {[round(c, 2) for c in ch[-4:]]}: {"converged" if done else "not yet"}')
                continue
            pid = alive(p + 'run.pid', inp)
            if not done or pid is None:
                continue
            os.kill(pid, signal.SIGTERM)
            for _ in range(600):
                if alive(p + 'run.pid', inp) is None:
                    break
                time.sleep(1)
            if alive(p + 'run.pid', inp) is None:
                with open(p + 'run.exit', 'w') as f:
                    f.write('0\n')
                note = (f'== {now()} converged (tools/converge_stop.py): the last recurrence intervals '
                        f'{[round(r, 2) for r in rec[-o.cycles - 1:]]} yr changed by less than {o.tol} yr')
                with open(p + 'run.log', 'a') as f:
                    f.write(note + '\n')
                print(f'{now()} {name}: {note[3:]}; marked finished', flush=True)
                left -= 1
        if o.dry_run:
            return
        if not widened and batch_finished(os.path.abspath(o.widen_after)):
            print(f'{now()} {o.widen_after} has finished: restarting {batch} with {o.jobs} jobs', flush=True)
            batch_sh = os.path.join(ROOT, 'tools', 'batch.sh')
            subprocess.run([batch_sh, 'stop', batch], cwd=ROOT)
            subprocess.run([batch_sh, 'start', batch, '-j', str(o.jobs)], cwd=ROOT)
            widened = True
        if left == 0 and widened:
            print(f'{now()} no watched run is left; done', flush=True)
            return
        time.sleep(o.interval)


if __name__ == '__main__':
    sys.exit(main())
