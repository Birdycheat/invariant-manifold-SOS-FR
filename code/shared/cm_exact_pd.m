function [ok, info] = cm_exact_pd(Q)
%CM_EXACT_PD Prove rational matrix positivity by exact fraction-free elimination.
Q = cm_exact(Q);
n = size(Q, 1);
info = struct('method', 'exact Bareiss leading principal minors', ...
    'size', n, 'positive_pivots', 0, 'failure_index', 0);
ok = false;
if n == 0 || size(Q, 2) ~= n || ~isequal(Q, Q.') || ~isempty(symvar(Q)), return; end
if ~all(isSymType(Q(:), 'integer') | isSymType(Q(:), 'rational')), return; end
[~, denominators] = numden(Q);
% Clear denominators with a positive factor before fraction-free elimination.
denominator = lcm(denominators(:));
A = Q * denominator;
previous = sym(1);
% Bareiss pivots are leading principal minors; strict positivity proves PD.
for k = 1:n
    pivot = A(k, k);
    if ~isAlways(pivot > 0)
        info.failure_index = k;
        return;
    end
    info.positive_pivots = k;
    if k < n
        remainder = k+1:n;
        A(remainder, remainder) = (A(remainder, remainder)*pivot ...
            - A(remainder, k)*A(k, remainder)) / previous;
    end
    previous = pivot;
end
ok = true;
end
