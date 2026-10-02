function out = h5read_reshape(fileName, datasetName, a)
% Read an SCycle dataset and give it a useful shape.
%
%   a.Ny, a.Nz       grid size (from domain.txt)
%   a.sI, a.eI       first and last saved step to read (time series only)
%   a.stride         read every stride-th saved step
%
% Context datasets (no time axis): body fields of size Ny*Nz become [Nz, Ny]
% (index Ii = iy*Nz + iz, so rows are depth and columns distance from the
% fault); other fields become column vectors.
% Time series (written with PETSc timestepping, HDF5 dims [Nt, N, 1]):
% body fields become [Nz, Ny, nt]; fault and boundary fields [N, nt];
% scalars per step become [nt, 1] (time, dt, regime, SS_index).

info = h5info(fileName, datasetName);
sz = info.Dataspace.Size;

if numel(sz) < 3 % context field
  data = h5read(fileName, datasetName);
  if numel(data) == a.Ny*a.Nz && a.Ny > 1 && a.Nz > 1
    out = reshape(data, a.Nz, a.Ny);
  else
    out = data(:);
  end
  return
end

Nt = sz(end);
eI = min(a.eI, Nt);
sI = max(1, min(a.sI, eI));
stride = max(1, a.stride);
nt = floor((eI - sI)/stride) + 1;

start = ones(1, numel(sz)); start(end) = sI;
count = sz; count(end) = nt;
strd = ones(1, numel(sz)); strd(end) = stride;
data = h5read(fileName, datasetName, start, count, strd);

N = prod(count(1:end-1));
if N == a.Ny*a.Nz && a.Ny > 1 && a.Nz > 1
  out = reshape(data, a.Nz, a.Ny, nt);
elseif N == 1 && (startsWith(datasetName,'/time/') || startsWith(datasetName,'/steadyState/'))
  out = reshape(data, nt, 1);
else
  out = reshape(data, N, nt);
end
end
