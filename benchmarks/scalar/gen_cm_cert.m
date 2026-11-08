function [ok, mu, info] = gen_cm_cert(P, p, r, ops)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
N = P.N; m = P.m;
phi = cell(N,1);
for i = 1:N, phi{i} = [0 0 P.H(i,2:p)]; end
xa = 0;
for i = 1:N, xa = padd(xa, P.a(i)*phi{i}); end
F = padd(P.g, pshift(xa, 1));
Ni = cell(N,1);
for i = 1:N
    Bphi = 0;
    for j = 1:N, Bphi = padd(Bphi, P.B(i,j)*phi{j}); end
    Ni{i} = padd(padd(conv(pder(phi{i}), F), -P.Lam(i,i)*phi{i}), ...
                 padd(-pshift(Bphi,1), [0 0 -P.v(i) -P.w(i)]));
end
% Factor the prescribed central order and the transverse invariance defect.
Fd = pdiv(F, m);
info.res_order = min(cellfun(@(c) find(abs(c) > 1e-9, 1) - 1, Ni));
x = sdpvar(1); mu = sdpvar(1);
toy = @(c) c*monolist(x, numel(c)-1);
% Completing squares gives a univariate matrix-positivity problem on the strip.
G = sdpvar(N+1, N+1);
G(1,1) = toy(-Fd);
for i = 1:N
    G(1,1+i) = -(P.a(i) - 2*mu*toy(pdiv(Ni{i}, m-1)))/2;
    G(1+i,1) = G(1,1+i);
    for j = i:N
        Mij = padd(padd(P.Lam(i,j)*(i==j), [0 P.B(i,j)]), -P.a(j)*pshift(pder(phi{i}),1));
        Mji = padd(padd(P.Lam(j,i)*(i==j), [0 P.B(j,i)]), -P.a(i)*pshift(pder(phi{j}),1));
        G(1+i,1+j) = -mu*(toy(Mij) + toy(Mji));
        G(1+j,1+i) = G(1+i,1+j);
    end
end
degx = max(arrayfun(@(k) degree(G(k), x), 1:numel(G)));
dS = max(0, 2*ceil(degx/2) - 2);
info.dS = dS;
S1 = sdpvar(N+1, N+1); cS = [];
for i = 1:N+1
    for j = i:N+1
        [s, c] = polynomial(x, dS);
        S1(i,j) = s; S1(j,i) = s; cS = [cS; c];
    end
end
Fc = [sos(G - 1e-3*eye(N+1) - (r^2 - x^2)*S1), sos(S1), mu >= 0, mu <= 1e4];
t = tic;
sol = solvesos(Fc, [], ops, [mu; cS]);
info.time = toc(t);
info.info = sol.info;
ok = (sol.problem == 0);
mu = value(mu);
end

function r = padd(p, q)
n = max(numel(p), numel(q));
r = [p zeros(1, n-numel(p))] + [q zeros(1, n-numel(q))];
end
function r = pshift(p, k)
r = [zeros(1,k) p];
end
function r = pder(p)
if numel(p) <= 1, r = 0; else, r = p(2:end).*(1:numel(p)-1); end
end
function r = pdiv(p, k)
scale = max(1, max(abs(p)));
assert(all(abs(p(1:min(k,numel(p)))) <= 1e-8*scale), 'division par x^%d non exacte', k);
r = p(k+1:end);
if isempty(r), r = 0; end
end
