function [ok, info] = cm_verify_exact(P, mu, p, r)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
x = sym('x');
N = P.N; m = P.m;
[Lam, B, a, v, w] = exact_data(P);
% Fix rational parameters, rounding the radius down before the proof.
mus = sym(round(1e6*mu))/1e6; rs = sym(floor(1e4*r))/1e4;
info.mu = mus; info.r = rs;
t0 = tic;
g0 = sum(P.gS.*x.^(0:numel(P.gS)-1));
phi = P.HS(:,2:p)*x.^((2:p).');
dphi = diff(phi, x);
F = expand(g0 + x*(a.'*phi));
[qF, rF] = quorem(F, x^m, x);
R = expand(dphi*F - Lam*phi - x*(B*phi) - v*x^2 - w*x^3);
qR = sym(zeros(N,1)); rR = sym(zeros(N,1));
for i = 1:N, [qR(i), rR(i)] = quorem(R(i), x^(m-1), x); end
info.divisible = isequal(rF, sym(0)) && isequal(rR, sym(zeros(N,1)));
ok = false; info.fail = 0;
if ~info.divisible, info.time = toc(t0); return; end
M = Lam + x*B - x*dphi*a.';
G = [-qF, -(a - 2*mus*qR).'/2; -(a - 2*mus*qR)/2, -mus*(M + M.')];
info.G = G;
info.deg = zeros(N+1, 1);
% Positive Sylvester minors cannot change sign without a root in the closed strip.
for k = 1:N+1
    Dk = expand(det(G(1:k,1:k)));
    info.deg(k) = polynomialDegree(Dk, x);
    if ~isAlways(subs(Dk, x, 0) > 0) || cm_count_roots(Dk, x, -rs, rs) > 0
        info.fail = k; info.time = toc(t0); return
    end
end
ok = true;
info.time = toc(t0);
end

