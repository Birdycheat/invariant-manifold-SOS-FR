function [ok, S, info] = jp_cert(P, d2, R, ops)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
N = P.N;
x = sdpvar(1); y = sdpvar(N,1);
f = [P.g*monolist(x, numel(P.g)-1) + x*P.a'*y; P.Lam*y + x*P.B*y + x^2*P.v + x^3*P.w];
mx = monolist(x, d2);
cJ = sdpvar(d2, 1);  J = cJ'*mx(2:end);
cH = sdpvar(N, d2);  H = cH*mx(1:d2);
Pq = sdpvar(N, N);
% Keep the specified partially quadratic dependence on the transverse variables.
V = J + y'*H + y'*Pq*y;
dV = jacobian(V, [x;y])*f;
[s, cs] = polynomial([x;y], 2*floor((degree(dV)-2)/2));
% This adaptation uses a ball and no strict decrease margin.
Fc = [sos(V - 1e-5*(x^2 + y'*y)), sos(-dV - s*(R^2 - x^2 - y'*y)), sos(s)];
prm = [cJ; cH(:); Pq(find(triu(ones(N)))); cs];
t = tic;
[sol, ~, ~, res] = solvesos(Fc, [], ops, prm);
info.time = toc(t); info.info = sol.info; info.res = max(abs(res)); info.sol = sol;
ok = (sol.problem == 0);
S.cJ = value(cJ); S.cH = value(cH); S.Pq = value(Pq);
[info.E, info.cval] = monomials(S, N);
end

function [E, c] = monomials(S, N)
d2 = numel(S.cJ);
E = zeros(0, N+1); c = zeros(0, 1);
for k = 1:d2
    E(end+1,:) = [k, zeros(1,N)]; c(end+1,1) = S.cJ(k);
end
for i = 1:N
    for k = 1:d2
        e = zeros(1,N); e(i) = 1;
        E(end+1,:) = [k-1, e]; c(end+1,1) = S.cH(i,k);
    end
end
for i = 1:N
    for j = i:N
        e = zeros(1,N); e(i) = e(i) + 1; e(j) = e(j) + 1;
        E(end+1,:) = [0, e]; c(end+1,1) = S.Pq(i,j)*(1 + (i ~= j));
    end
end
end
