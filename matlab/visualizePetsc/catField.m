function out = catField(out, fieldName, dim, data)
% Append data to out.(fieldName) along dimension dim, or create the field.
% Useful to join output of a run that was restarted into a new outputDir.
if isfield(out, fieldName) && ~isempty(out.(fieldName))
  out.(fieldName) = cat(dim, out.(fieldName), data);
else
  out.(fieldName) = data;
end
end
