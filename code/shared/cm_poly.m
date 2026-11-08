function value = cm_poly(operation, varargin)
%CM_POLY Sparse multivariate polynomials with exact rational coefficients.
switch lower(operation)
    case 'zero'
        value = struct('E', zeros(0, varargin{1}), 'c', sym(zeros(0, 1)));
    case 'mono'
        coefficient = sym(1);
        if numel(varargin) > 1, coefficient = cm_exact(varargin{2}); end
        value = normalise(double(varargin{1}(:).'), coefficient);
    case 'from'
        value = normalise(double(varargin{1}), cm_exact(varargin{2}));
    case 'add'
        a = varargin{1}; b = varargin{2}; scale = sym(1);
        if numel(varargin) > 2, scale = cm_exact(varargin{3}); end
        value = normalise([a.E; b.E], [a.c; scale*b.c]);
    case 'scale'
        a = varargin{1}; value = normalise(a.E, a.c*cm_exact(varargin{2}));
    case 'mul'
        a = varargin{1}; b = varargin{2}; nv = size(a.E, 2);
        if isempty(a.c) || isempty(b.c), value = cm_poly('zero', nv); return; end
        exponents = reshape(reshape(a.E, [], 1, nv)+reshape(b.E, 1, [], nv), [], nv);
        coefficients = a.c*b.c.';
        value = normalise(exponents, coefficients(:));
    case 'deriv'
        a = varargin{1}; j = varargin{2}; keep = a.E(:, j) > 0;
        if ~any(keep), value = cm_poly('zero',size(a.E,2)); return; end
        exponents = a.E(keep, :); coefficients = a.c(keep).*exponents(:, j);
        exponents(:, j) = exponents(:, j)-1;
        value = normalise(exponents, coefficients);
    case 'eval'
        a = varargin{1}; point = cm_exact(varargin{2});
        if isempty(a.c), value = sym(0); return; end
        monomials = sym(ones(size(a.c)));
        for j = 1:size(a.E, 2), monomials = monomials.*point(j).^a.E(:, j); end
        value = sum(a.c.*monomials);
    % Collect all monomial products in v'Qv, including both off-diagonal terms.
    case 'gram'
        E = double(varargin{1}); Q = cm_exact(varargin{2}); nv = size(E, 2);
        assert(isequal(size(Q), [size(E,1), size(E,1)]));
        exponents = reshape(reshape(E, [], 1, nv)+reshape(E, 1, [], nv), [], nv);
        value = normalise(exponents, Q(:));
    case 'equal'
        difference = cm_poly('add', varargin{1}, varargin{2}, -1);
        value = isempty(difference.c);
    case 'compose'
        a = varargin{1}; replacement = varargin{2}; nv = size(replacement{1}.E, 2);
        value = cm_poly('zero', nv);
        for k = 1:numel(a.c)
            term = cm_poly('mono', zeros(1,nv), a.c(k));
            for j = 1:size(a.E, 2)
                for power = 1:a.E(k,j), term = cm_poly('mul', term, replacement{j}); end
            end
            value = cm_poly('add', value, term);
        end
    case 'truncate'
        a = varargin{1}; keep = sum(a.E,2) <= varargin{2};
        value = normalise(a.E(keep,:), a.c(keep));
    case 'expression'
        a = varargin{1}; variables = varargin{2};
        if isempty(a.c), value = sym(0); return; end
        terms = a.c;
        for j = 1:size(a.E,2), terms = terms.*variables(j).^a.E(:,j); end
        value = sum(terms);
    otherwise
        error('cm:poly:operation', 'Unknown polynomial operation: %s', operation);
end
end

function value = normalise(E, coefficients)
coefficients = coefficients(:);
assert(size(E,1) == numel(coefficients) && all(E(:) >= 0) && all(E(:) == fix(E(:))));
if isempty(coefficients)
    value = struct('E', E, 'c', sym(zeros(0,1)));
    return;
end
% Combine identical exponents exactly so a zero polynomial has no terms.
[exponents, ~, groups] = unique(E, 'rows');
[~,order] = sort(groups);
endpoints = cumsum(accumarray(groups,1,[size(exponents,1),1]));
partialSums = cumsum(coefficients(order));
combined = partialSums(endpoints);
if numel(combined)>1, combined(2:end) = combined(2:end)-partialSums(endpoints(1:end-1)); end
keep = ~isAlways(combined == 0);
value = struct('E', exponents(keep,:), 'c', combined(keep));
end
