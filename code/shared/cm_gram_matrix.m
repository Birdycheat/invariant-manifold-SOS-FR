function matrix = cm_gram_matrix(value, form)
%CM_GRAM_MATRIX Decode a full or row-wise upper-triangular exact Gram matrix.
if nargin < 2, form = 'triangular'; end
if strcmp(form, 'full'), matrix = cm_exact(value); return; end
n = numel(value);
matrix = sym(zeros(n));
for k = 1:n
    if iscell(value), row = cm_exact(value{k}); else, row = cm_exact(value(k,:)); end
    row = row(:).';
    assert(numel(row) == n-k+1, 'cm:gram:shape', 'Invalid triangular Gram row.');
    matrix(k,k:n) = row;
    matrix(k:n,k) = row.';
end
end
