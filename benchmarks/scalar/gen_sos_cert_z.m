function [ok, info] = gen_sos_cert_z(P, D, r, ops, p, epsm)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if nargin < 6, epsm = 1e-3; end
x = sdpvar(1); z = sdpvar(P.N,1);
mx = monolist(x, p);
if p >= 2
    phi = P.H(:,2:p)*mx(3:p+1);
    dphi = P.H(:,2:p)*((2:p)'.*mx(2:p));
else
    phi = zeros(P.N, 1); dphi = zeros(P.N, 1);
end
g0 = P.g*monolist(x, numel(P.g)-1);
xd = g0 + x*P.a'*(phi + z);
% Apply the chain rule for z = y - phi(x), including -Dphi*xd.
f = [xd; P.Lam*(phi + z) + x*P.B*(phi + z) + x^2*P.v + x^3*P.w - dphi*xd];
% Fix the degree in adapted coordinates before the with/without-FR comparison.
mons = monolist([x;z], D);
dg = zeros(length(mons),1);
for i = 1:length(mons), dg(i) = degree(mons(i)); end
mons = mons(dg >= 2);
c = sdpvar(length(mons),1);
W = c'*mons;
dW = jacobian(W, [x;z])*f;
k = min(2*P.m-2, D+P.m-1);
ds = 2*floor((degree(dW)-2)/2);
[s, cs] = polynomial([x;z], ds);
Fc = [sos(W), sos(-dW - epsm*x^k - s*(r^2 - x^2)), sos(s)];
t = tic;
[sol, vb, Qb, res] = solvesos(Fc, [], ops, [c; cs]);
info.time = toc(t); info.info = sol.info; info.res = max(abs(res));
ok = (sol.problem == 0);
info.v = vb; info.Q = Qb; info.x = x; info.z = z; info.k = k;
Wn = value(c)'*mons;
sn = replace(s, cs, value(cs));
info.rho = -jacobian(Wn, [x;z])*f - epsm*x^k - sn*(r^2 - x^2) - vb{2}'*Qb{2}*vb{2};
vars = [x; z]; E = zeros(length(mons), length(vars));
for i = 1:length(mons)
    for j = 1:length(vars), E(i,j) = degree(mons(i), vars(j)); end
end
info.E = E; info.cval = value(c);
end
