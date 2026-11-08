function [ok, out] = jet_sos_center(P, D, r, p, fn)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
ops = sdpsettings('solver','mosek','verbose',0);
N = P.N; m = P.m;
% Round the radius down before constructing the prescribed jet face.
rn = floor(1e4*r); r = rn/1e4;
[~, info0] = gen_sos_cert_z(P, D, r, ops, p);
x = info0.x; z = info0.z; vars = [x; z]; k = info0.k;
E = cell(1, 3);
for b = 1:3, E{b} = monexp(info0.v{b}, vars); end
low = @(Eb) Eb(:,1) <= (m-1)/2 & all(Eb(:,2:end) == 0, 2);
zdeg = @(Eb) sum(Eb(:,2:end), 2);
E{2} = E{2}(~low(E{2}) & zdeg(E{2}) <= D/2, :);
E{3} = E{3}(~low(E{3}) & zdeg(E{3}) <= D/2, :);
mx = monolist(x, p);
if p >= 2
    phi = P.H(:,2:p)*mx(3:p+1);
    dphi = P.H(:,2:p)*((2:p)'.*mx(2:p));
else
    phi = zeros(N, 1); dphi = zeros(N, 1);
end
g0 = P.g*monolist(x, numel(P.g)-1);
xd = g0 + x*P.a'*(phi + z);
f = [xd; P.Lam*(phi + z) + x*P.B*(phi + z) + x^2*P.v + x^3*P.w - dphi*xd];
out.deleted = {};
maxit = 10;
% Optional tiny-diagonal deletions retry numerical centering on smaller bases.
for it = 1:maxit
    [ok, out1] = solve_center(E, vars, f, x, k, r, ops);
    if ~ok || out1.t > 1e-7 || it == maxit, break; end
    nd = 0;
    for b = 1:3
        dg = diag(out1.Q{b});
        del = dg < 1e-8*max(dg);
        out.deleted = [out.deleted; arrayfun(@(i) sprintf('Q%d:%s', b, mat2str(E{b}(i,:))), ...
                       find(del)', 'uni', 0)'];
        E{b} = E{b}(~del, :); nd = nd + nnz(del);
    end
    if nd == 0, break; end
end
out = mergestruct(out, out1);
out.iter = it; out.E = E; out.r = r; out.k = k;
if nargin < 5 || isempty(fn), return; end
% Export ordered bases and numerical Gram blocks for native rational reconstruction.
[Lam, B, a, vv, w] = exact_data(P);
S = struct('N', N, 'm', m, 'p', p, 'D', D, 'k', k, 'r', sprintf('%d/10000', rn), ...
           'solver_ok', ok, 't', out.t, 'e', out.e, ...
           'Lam', rat(Lam), 'B', rat(B), 'a', rat(a), 'v', rat(vv), 'w', rat(w), ...
           'g', rat(P.gS), 'H', rat(P.HS(:,1:p)), ...
           'E1', flat(E{1}), 'E2', flat(E{2}), 'E3', flat(E{3}), ...
           'Q1', flat(out.Q{1}), 'Q2', flat(out.Q{2}), 'Q3', flat(out.Q{3}));
fid = fopen(fn, 'w'); fprintf(fid, '%s', jsonencode(S)); fclose(fid);
end

function [ok, out] = solve_center(E, vars, f, x, k, r, ops)
n = cellfun(@(Eb) size(Eb, 1), E);
v = cellfun(@(Eb) monvec(vars, Eb), E, 'uni', 0);
Q = arrayfun(@(nb) sdpvar(nb), n, 'uni', 0);
e = sdpvar(1); t = sdpvar(1);
W = v{1}'*Q{1}*v{1};
s = v{3}'*Q{3}*v{3};
T = -jacobian(W, vars)*f - e*x^k - s*(r^2 - x^2) - v{2}'*Q{2}*v{2};
% Fix the Gram scale before maximizing the common eigenvalue margin.
F = [coefficients(T, vars) == 0, e >= t, trace(Q{1}) + trace(Q{2}) + trace(Q{3}) == 1];
for b = 1:3, F = [F, Q{b} - t*eye(n(b)) >= 0]; end
tic; sol = optimize(F, -t, ops); out.time = toc;
ok = sol.problem == 0;
out.t = value(t); out.e = value(e); out.n = n;
out.Q = cellfun(@value, Q, 'uni', 0);
out.info = sol.info;
end

function a = mergestruct(a, b)
for fn = fieldnames(b)', a.(fn{1}) = b.(fn{1}); end
end

function Eb = monexp(vb, vars)
Eb = zeros(length(vb), length(vars));
for i = 1:length(vb)
    for j = 1:length(vars)
        try, Eb(i,j) = degree(vb(i), vars(j)); catch, end
    end
end
end

function v = monvec(vars, Eb)
v = [];
for i = 1:size(Eb, 1)
    mi = 1;
    for j = 1:length(vars), if Eb(i,j) > 0, mi = mi*vars(j)^Eb(i,j); end, end
    v = [v; mi];
end
end

function s = rat(A)
s = struct('dims', size(A), 'vals', {arrayfun(@(q) char(q), A(:)', 'uni', 0)});
end

function s = flat(A)
s = struct('dims', size(A), 'vals', {arrayfun(@(q) sprintf('%.17g', q), A(:)', 'uni', 0)});
end
