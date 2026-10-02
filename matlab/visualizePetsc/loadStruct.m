function dataStruct = loadStruct(fileName, delim)
% Loads the contents of the file (input: fileName) into a struct (output: dataStruct).
%
% The file is assumed to consist of 2 columns separated by the specified delimiter.
% For example, if delim = ' = ', then the file is assumed to have the structure:
%
% var1 = 1.0
% var2 = [0 15 30]
% var3 = RK43 # comments are formatted like this
%
% The data struct will contain fields with the names {var1,var2,...}. Numbers
% and [lists] become numeric values, anything else a character array.
% Lines without the delimiter are skipped; text after '#' is ignored.

dataStruct = struct();
fid = fopen(fileName);
if fid < 0
  error('loadStruct: cannot open %s', fileName);
end

while ~feof(fid)
  fileLine = fgetl(fid); % load current line
  if ~ischar(fileLine) || isempty(strtrim(fileLine)), continue, end

  % split line at the first delimiter
  matches = strfind(fileLine,delim);
  if isempty(matches), continue, end
  fieldName = strtrim(fileLine(1:matches(1)-1));
  fieldValue = fileLine(matches(1)+length(delim):end);

  % remove any trailing comments
  commentIndex = strfind(fieldValue,'#');
  if ~isempty(commentIndex)
    fieldValue = fieldValue(1:commentIndex(1)-1);
  end
  fieldValue = strtrim(fieldValue);

  value = str2double(fieldValue); % numbers
  if isnan(value) && startsWith(fieldValue,'[') && endsWith(fieldValue,']') % lists such as [0 15 30]
    value = sscanf(fieldValue(2:end-1),'%f')';
  end
  if isempty(value) || (isscalar(value) && isnan(value) && ~strcmpi(fieldValue,'nan'))
    value = fieldValue; % text
  end
  dataStruct.(matlab.lang.makeValidName(fieldName)) = value;
end
fclose(fid);

end
