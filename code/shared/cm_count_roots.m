function count = cm_count_roots(polynomial, variable, lower, upper)
%CM_COUNT_ROOTS Count distinct roots on a closed interval by exact Sturm signs.
polynomial = expand(polynomial);
assert(~isequal(polynomial, sym(0)), 'cm:roots:zero', 'The zero polynomial has infinitely many roots.');
assert(isAlways(lower <= upper), 'cm:roots:interval', 'Invalid interval.');
if polynomialDegree(polynomial, variable) == 0, count = 0; return; end
% Remove repeated factors so Sturm counts distinct roots.
[polynomial, remainder] = quorem(polynomial, gcd(polynomial, diff(polynomial, variable)), variable);
assert(isequal(remainder, sym(0)));
sequence = {polynomial, diff(polynomial, variable)};
while polynomialDegree(sequence{end}, variable) > 0
    [~, remainder] = quorem(sequence{end-1}, sequence{end}, variable);
    if isequal(remainder, sym(0)), break; end
    sequence{end+1} = -expand(remainder); %#ok<AGROW>
end
% Sturm variations include the upper endpoint; add a root at the lower one.
count = variations(sequence, variable, lower) - variations(sequence, variable, upper);
if isfinite(lower) && isAlways(subs(polynomial, variable, lower) == 0), count = count+1; end
end

function changes = variations(sequence, variable, point)
signs = zeros(1, numel(sequence));
for k = 1:numel(sequence)
    if isfinite(point)
        value = subs(sequence{k}, variable, point);
    else
        degree = double(polynomialDegree(sequence{k}, variable));
        coefficients = coeffs(sequence{k}, variable, 'All');
        value = coefficients(1);
        if point < 0, value = value*(-1)^degree; end
    end
    if isAlways(value > 0), signs(k) = 1;
    elseif isAlways(value < 0), signs(k) = -1;
    else, assert(isAlways(value == 0), 'cm:roots:sign', 'Unresolved exact sign.');
    end
end
signs = signs(signs ~= 0);
changes = nnz(signs(1:end-1) ~= signs(2:end));
end
