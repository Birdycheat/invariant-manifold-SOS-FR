function values = cm_exact_flat(data, binary)
%CM_EXACT_FLAT Decode a column-major dims/vals coefficient array.
if nargin < 2, binary = false; end
dimensions = double(data.dims(:).');
if isscalar(dimensions), dimensions = [dimensions, 1]; end
values = reshape(cm_exact(data.vals, binary), dimensions);
end
