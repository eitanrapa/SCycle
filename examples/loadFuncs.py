"""Load SCycle output (HDF5 files and the .txt context files) into nested dictionaries.

SCycle writes <outputDir>data_context.h5, data_1D.h5, data_2D.h5 and data_steadyState.h5,
plus <outputDir>domain.txt, mediator.txt, fault.txt, ... where outputDir is the prefix
given in the input file (for example "data/ex2_"). Pass either that prefix or the older
style without the trailing underscore ("data/ex2"); both work.

Each HDF5 path becomes nested keys: '/fault/slipVel' -> d['fault']['slipVel'].
By default every dataset in a file is loaded. With Ny and Nz given, fields are reshaped:
body fields to (Ny, Nz[, Nt]), fault and boundary fields to (Nz or Ny[, Nt]), time series
to (Nt,), with time as the last axis (index Ii = iy*Nz + iz, z fastest).
"""
import h5py
import numpy as np
import os.path


def _prefix(filePath):
  """Return the output prefix that SCycle prepended to its file names."""
  if filePath.endswith('_') or filePath.endswith('/') or filePath.endswith(os.sep):
    return filePath
  if os.path.isfile(filePath + 'data_context.h5') or os.path.isfile(filePath + 'domain.txt'):
    return filePath
  return filePath + '_'


def _datasetKeys(currFile):
  keys = []
  currFile.visititems(lambda name, obj: keys.append('/' + name) if isinstance(obj, h5py.Dataset) else None)
  return keys


def _loadFile(fileName, keys, d, Ny=-1, Nz=-1):
  with h5py.File(fileName, 'r') as currFile:
    if isinstance(keys, str) and keys == 'default':
      keys = _datasetKeys(currFile)
    for key in keys:
      vec = loadDataSet(currFile, key)
      if vec.size != 0:
        vec = reshapeDataSet(vec, Ny, Nz)
        createSubDicts(d, key, vec)
  return d


def loadContext(filePath, keys='default'):
  """Mesh, material and fault parameters (data_context.h5) plus the .txt context files."""
  prefix = _prefix(filePath)
  d = {}
  for name in ['domain', 'mediator', 'heatEquation', 'fault', 'momBal', 'pressureEq', 'grainSizeEv']:
    txtFileName = "%s%s.txt" % (prefix, name)
    if os.path.isfile(txtFileName):
      appendDict(d, 'med' if name == 'mediator' else name, load_txt(txtFileName))

  Ny = d.get('domain', {}).get('Ny', -1)
  Nz = d.get('domain', {}).get('Nz', -1)
  return _loadFile('%sdata_context.h5' % prefix, keys, d, Ny, Nz)


def load_1D(filePath, indict=None, keys='default', Nz=-1, Ny=-1):
  """Fault and boundary time series (data_1D.h5): fields of size Nz or Ny per saved step."""
  d = {} if indict is None else indict
  return _loadFile('%sdata_1D.h5' % _prefix(filePath), keys, d, Ny, Nz)


def load_2D(filePath, indict=None, keys='default', Nz=-1, Ny=-1):
  """Body fields (data_2D.h5): fields of size Ny*Nz per saved step."""
  d = {} if indict is None else indict
  return _loadFile('%sdata_2D.h5' % _prefix(filePath), keys, d, Ny, Nz)


def load_SS(filePath, indict=None, keys='default', Nz=-1, Ny=-1):
  """Steady-state iterations (data_steadyState.h5)."""
  d = {} if indict is None else indict
  return _loadFile('%sdata_steadyState.h5' % _prefix(filePath), keys, d, Ny, Nz)


def load_txt(txtFileName):
  """Parse 'key = value # comment' lines into a dict, converting numbers and [lists]."""
  d = {}
  with open(txtFileName, 'r') as currFile:
    for currLine in currFile.readlines():
      currStrList = currLine.split(' = ', 1)
      if len(currStrList) == 2:
        key = currStrList[0].strip()
        val = currStrList[1].split('#')[0].strip()  # drop the comment and the newline
        str2Value(d, key, val)
  return d


def createSubDicts(d, key, vec):
  """Create nested dictionaries for the parts of an HDF5 path and store vec at the leaf."""
  currd = d
  keylist = key.split('/')
  for subkey in keylist[:-1]:
    if subkey != '':
      if subkey not in currd.keys():
        currd[subkey] = {}
      currd = currd[subkey]
  currd[keylist[-1]] = vec


def appendDict(d, key, d2):
  if key in d.keys():
    d[key] = dict(d[key], **d2)
  else:
    d[key] = d2


def loadDataSet(currFile, key):
  try:
    return np.array(currFile[key])
  except KeyError:
    return np.array([])


def reshapeDataSet(vec, Ny, Nz):
  """Reshape a flat or (Nt, N, 1) dataset to (Ny, Nz[, Nt]) or (N[, Nt]) with time last."""
  if not (Nz > -1 and Ny > -1 and isinstance(vec, np.ndarray)):
    return vec
  if vec.ndim == 1:  # context fields, no time axis
    if vec.size == Ny * Nz and Ny > 1 and Nz > 1:
      return vec.reshape(Ny, Nz)
    return vec
  Nt = vec.shape[0]
  N = vec.size // Nt if Nt > 0 else 0
  if N == Ny * Nz and Ny > 1 and Nz > 1:
    return np.moveaxis(vec.reshape(Nt, Ny, Nz), 0, -1)
  if N == 1:
    return vec.reshape(Nt)
  return np.moveaxis(vec.reshape(Nt, N), 0, -1).squeeze()


def str2FloatList(insStr):
  """Convert '[a b c]' to a list of floats; return False if insStr is not such a list."""
  if '[' in insStr and ']' in insStr:
    strList = insStr.strip().lstrip('[').rstrip(']').split()
    try:
      return [float(x) for x in strList]
    except ValueError:
      return False
  return False


def str2Value(d, key, val):
  temp = str2FloatList(val)
  if not isinstance(temp, bool):
    d[key] = temp
  else:
    try:
      d[key] = int(val)
    except ValueError:
      try:
        d[key] = float(val)
      except ValueError:
        d[key] = val
