function P = gen_instance_exact(N, m, d, seed, K, wscale, ascale)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if nargin < 5, K = []; end
if nargin < 6, wscale = []; end
if nargin < 7, ascale = 1; end
P = gen_instance(N, m, d, seed, K, wscale, ascale);
% Freeze the sampled coefficients on a rational grid before recomputing the jet.
D3 = @(c) sym(round(1000*c))/1000;
P.LamS = D3(P.Lam); P.BS = D3(P.B); P.aS = D3(P.a); P.vS = D3(P.v); P.wS = D3(P.w);
P.Lam = double(P.LamS); P.lam = -diag(P.Lam); P.B = double(P.BS);
P.a = double(P.aS); P.v = double(P.vS); P.w = double(P.wS);
% Recompute criticality and invariance using exact rational arithmetic.
[P.gS, P.HS] = exact_jet(P, size(P.H, 2), true);
P.g = double(P.gS);
P.H = double(P.HS);
end
