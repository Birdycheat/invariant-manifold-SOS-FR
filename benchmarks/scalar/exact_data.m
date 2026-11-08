function [Lam, B, a, v, w] = exact_data(P)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
S = @(c) sym(c, 'f');
% Prefer retained rational coefficients; otherwise interpret stored doubles exactly.
if isfield(P, 'LamS')
    Lam = P.LamS; B = P.BS; a = P.aS; v = P.vS; w = P.wS;
else
    Lam = S(P.Lam); B = S(P.B); a = S(P.a); v = S(P.v); w = S(P.w);
end
end
