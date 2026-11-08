function rows = cm_json_rows(value)
%CM_JSON_ROWS Normalise homogeneous or mixed JSON rows to a cell matrix.
if isnumeric(value)
    rows = num2cell(value);
elseif iscell(value) && isvector(value) && ~isempty(value) && iscell(value{1})
    width = numel(value{1});
    rows = cell(numel(value), width);
    for k = 1:numel(value), rows(k,:) = reshape(value{k}, 1, []); end
elseif iscell(value)
    rows = value;
else
    rows = num2cell(value);
end
end
