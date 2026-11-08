function q = cm_exact(value, binary)
%CM_EXACT Convert stored coefficients to exact symbolic rationals.
if nargin < 2, binary = false; end
if isa(value, 'sym'), q = value; return; end
if iscell(value)
    if isvector(value) && ~isempty(value) && iscell(value{1}), value = cm_json_rows(value); end
    if isempty(value), q = sym(zeros(size(value))); return; end
    if all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)),value(:)))
        q = parseArray(string(value));
    else
        q = sym(zeros(size(value)));
        for k = 1:numel(value), q(k) = cm_exact(value{k}, binary); end
    end
elseif isnumeric(value) || islogical(value)
    assert(isreal(value) && all(isfinite(value(:))), 'cm:exact:finite', ...
        'Coefficients must be finite and real.');
    % Preserve binary64 values exactly when checking floating candidates.
    if binary
        q = sym(double(value), 'f');
    else
        q = parseArray(compose('%.17g',double(value)));
    end
else
    q = parseArray(string(value));
end
end

function q = parseArray(texts)
if isempty(texts), q = sym(zeros(size(texts))); return; end
fragments = cell(numel(texts),1);
for k = 1:numel(texts), fragments{k} = rationalText(char(texts(k))); end
q = reshape(sym(fragments),size(texts));
assert(all(isAlways(isfinite(q(:)))), 'cm:exact:finite', 'Invalid rational coefficient.');
end

function text = rationalText(text)
text = strtrim(text);
assert(~isempty(regexp(text, '^[+-]?\d+(/\d+)?$', 'once')) || ...
    ~isempty(regexp(text, '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$', 'once')), ...
    'cm:exact:rational', 'Expected an integer, fraction, or decimal coefficient.');
if contains(text, '/')
    return;
elseif contains(lower(text), 'e') || contains(text, '.')
    pieces = split(string(lower(text)), 'e');
    mantissa = char(pieces(1));
    direction = '';
    if mantissa(1) == '-', direction = '-'; mantissa = mantissa(2:end);
    elseif mantissa(1) == '+', mantissa = mantissa(2:end); end
    position = strfind(mantissa, '.');
    decimals = 0;
    if ~isempty(position), decimals = numel(mantissa) - position; end
    digitsText = erase(mantissa, '.');
    exponent = 0;
    if numel(pieces) > 1, exponent = str2double(pieces(2)); end
    scale = exponent-decimals;
    if scale >= 0
        text = [direction,digitsText,repmat('0',1,scale)];
    else
        text = [direction,digitsText,'/1',repmat('0',1,-scale)];
    end
end
end
