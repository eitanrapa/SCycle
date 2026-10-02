function [groupNames, datasetNames] = recursiveGroupFinder(varargin)
% Full paths of all groups and datasets in an HDF5 file, from h5info output.
%
%   [groupNames, datasetNames] = recursiveGroupFinder(info, info.Groups, {})
%   [groupNames, datasetNames] = recursiveGroupFinder(info.Groups)
%
% groupNames and datasetNames are cell arrays of paths such as '/fault' and
% '/fault/slipVel'. The three-argument form also lists datasets at the root.

if nargin == 1
  info = []; groups = varargin{1};
else
  info = varargin{1}; groups = varargin{2};
end

groupNames = {};
datasetNames = {};
if ~isempty(info) && isfield(info,'Datasets') && ~isempty(info.Datasets)
  for d = 1:numel(info.Datasets)
    datasetNames{end+1} = ['/' info.Datasets(d).Name]; %#ok<AGROW>
  end
end

for g = 1:numel(groups)
  groupNames{end+1} = groups(g).Name; %#ok<AGROW>
  for d = 1:numel(groups(g).Datasets)
    datasetNames{end+1} = [groups(g).Name '/' groups(g).Datasets(d).Name]; %#ok<AGROW>
  end
  [subGroups, subDatasets] = recursiveGroupFinder(groups(g).Groups);
  groupNames = [groupNames subGroups]; %#ok<AGROW>
  datasetNames = [datasetNames subDatasets]; %#ok<AGROW>
end
end
