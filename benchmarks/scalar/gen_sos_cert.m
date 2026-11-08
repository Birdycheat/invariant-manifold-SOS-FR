function [ok, Wn, info] = gen_sos_cert(P, D, r, ops)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
x = sdpvar(1); y = sdpvar(P.N,1);
g0 = P.g*monolist(x, numel(P.g)-1);
f = [g0 + x*P.a'*y; P.Lam*y + x*P.B*y + x^2*P.v + x^3*P.w];
mons = monolist([x;y], D);
dg = zeros(length(mons),1);
for i = 1:length(mons), dg(i) = degree(mons(i)); end
% Exclude constant and linear terms from the Lyapunov candidate.
mons = mons(dg >= 2);
c = sdpvar(length(mons),1);
W = c'*mons;
dW = jacobian(W, [x;y])*f;
k = min(2*P.m-2, D+P.m-1);
ds = 2*floor((degree(dW)-2)/2);
[s, cs] = polynomial([x;y], ds);
% The localizer bounds only x; transverse coordinates remain unbounded.
Fc = [sos(W), sos(-dW - 1e-3*x^k - s*(r^2 - x^2)), sos(s)];
t = tic;
[sol, vb, Qb, res] = solvesos(Fc, [], ops, [c; cs]);
info.v = vb; info.Q = Qb;
info.time = toc(t);
info.info = sol.info;
info.res = max(abs(res));
info.sol = sol;
info.nmon = length(mons);
ok = (sol.problem == 0);
Wn = clean(value(c)'*mons, 1e-9);
info.x = x; info.y = y; info.f = f;
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%d_%d', P.N, D);
if ~isKey(cache, key)
    vars = [x; y]; E = zeros(length(mons), length(vars));
    for i = 1:length(mons)
        for j = 1:length(vars), E(i,j) = degree(mons(i), vars(j)); end
    end
    cache(key) = E;
end
info.E = cache(key);
info.cval = value(c);
end
