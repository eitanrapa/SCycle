#!/usr/bin/env python3
"""Flag keys in SCycle input files that no part of the code reads.

SCycle silently ignores unknown keys, so a typo (retartFromChkptSS) or a key
retired since an input was written (problemType) changes nothing and gives no
warning. This script collects every key the C++ sources compare
against (var.compare("key") / var == "key") and lists the keys of each input
file that are not among them.

Limits: a key is "known" if any component reads it, so a key meant for one
component but read only by another is not flagged. Example: linSolver is read
only by the pressure equation; the momentum balance reads linSolverSS and
linSolverTrans, and the heat equation linSolver_heateq. Check the component's loadSettings when in doubt.

Interior faults: for each name in interiorFaults = [...], the keys <name>_y and
<name>_<key> (for any key above) are accepted. Bulk state fields (hard_, ...): their common
keys, built at run time from the prefix, are added from each field's constructor.

Usage: tools/checkkeys.py file.in [file.in ...]   (exit status 1 if any are unknown)
"""
import os, re, sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(root, 'source')
known = set()
pat = re.compile(r'var\s*(?:\.compare\(\s*"([^"]+)"\s*\)|==\s*"([^"]+)")')
for name in os.listdir(src):
    if name.endswith('.cpp'):
        with open(os.path.join(src, name), errors='replace') as f:
            for line in f:
                if line.lstrip().startswith(('//', '//~')):
                    continue
                for a, b in pat.findall(line):
                    known.add(a or b)
# keys built at run time: DislocationCreep reads "disl" + prefix + "_AVals" for prefix "" and "2"
for m in re.finditer(r'"disl"\s*\+\s*_prefix\s*\+\s*"(_\w+)"', open(os.path.join(src, 'dislocationCreep.cpp')).read()):
    for prefix in ('', '2'):
        known.add('disl' + prefix + m.group(1))

# bulk state fields (bulkStateField.cpp): BulkStateField(D,"name","prefix","symbol") reads the common
# keys <prefix>type, <prefix><symbol>Vals, ... that parseCommon builds at run time
common = open(os.path.join(src, 'bulkStateField.cpp')).read()
plain = re.findall(r'_prefix\s*\+\s*"(\w+)"', common)
withsym = re.findall(r'_prefix\s*\+\s*_symbol\s*\+\s*"(\w+)"', common)
for name in os.listdir(src):
    if name.endswith('.cpp'):
        for m in re.finditer(r'BulkStateField\(D,\s*"(\w+)",\s*"(\w+)",\s*"(\w+)"\)', open(os.path.join(src, name), errors='replace').read()):
            prefix, symbol = m.group(2), m.group(3)
            known.update(prefix + k for k in plain)
            known.update(prefix + symbol + k for k in withsym)


def interior_faults(path):
    """Names listed by interiorFaults = [name ...] in an input file."""
    names = []
    with open(path) as f:
        for line in f:
            line = line.split('#')[0]
            if line.startswith('interiorFaults = '):
                names = line.split(' = ', 1)[1].strip().lstrip('[').rstrip(']').split()
    return names

status = 0
for path in sys.argv[1:]:
    unknown = []
    faults = interior_faults(path)
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.split('#')[0]
            if ' = ' not in line or line[0].isspace():
                continue
            key = line.split(' = ')[0].strip()
            # an interior fault <name> reads <name>_y and <name>_<key> for any key a fault reads
            # (the check accepts any known key there)
            for name in faults:
                if key == name + '_y' or (key.startswith(name + '_') and key[len(name) + 1:] in known):
                    key = None
                    break
            if key and key not in known:
                unknown.append((n, key))
    if unknown:
        status = 1
        print('%s: %d key(s) not read by SCycle' % (path, len(unknown)))
        for n, key in unknown:
            print('  line %d: %s' % (n, key))
    else:
        print('%s: all keys are read by SCycle' % path)
sys.exit(status)
