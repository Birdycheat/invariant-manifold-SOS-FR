function out = coupled_sos(action, varargin)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
out = feval(['cs_' action], varargin{:});
end

function P = cs_data(S)
N = S.num; P.ny = S.ny; P.r = N.r;
P.Fhat = reshape(N.Fhat, 2, 6); P.Ehat = reshape(N.Ehat, S.ny, 5);
P.a = reshape(N.a, 2, S.ny); P.Lam = reshape(N.Lam, S.ny, S.ny);
P.B1 = reshape(N.B1, S.ny, S.ny); P.B2 = reshape(N.B2, S.ny, S.ny);
P.phi = reshape(N.phi, S.ny, 3); P.eta4 = reshape(N.eta4, S.ny, 5);
end

function [vars, f, x] = cs_system(P, coords, enc)
x = sdpvar(2, 1); w = sdpvar(P.ny, 1); vars = [x; w];
m4 =[x(1)^4; x(1)^3*x(2); x(1)^2*x(2)^2; x(1)*x(2)^3; x(2)^4];
m5 = [x(1)^5; x(1)^4*x(2); x(1)^3*x(2)^2; x(1)^2*x(2)^3; x(1)*x(2)^4; x(2)^5];
phi = P.phi(:,1)*x(1)^2 + P.phi(:,2)*x(1)*x(2) + P.phi(:,3)*x(2)^2;
switch coords
    case 'xy',  z = w - phi;
    case 'xz',  z = w;
    case 'xz4', z = w + P.eta4*m4;
end
% Bound cosine and exponential factors by four polynomial enclosure vertices.
u = x(1) - x(2); rho = x'*x;
if enc == 1
    G1 = {1, 1 + u^2/2}; G2 = {1 - rho, 1};
else
    G1 = {1 + u^2/2 - u^4/24, 1 + u^2/2}; G2 = {1 - rho, 1 - rho + rho^2/2};
end
Fh = P.Fhat*m5; Eh = P.Ehat*m4;
xb = [x(1)*(P.a(1,:)*z); x(2)*(P.a(2,:)*z)];
zb = (P.Lam + x(1)*P.B1 + x(2)*P.B2)*z;
f = cell(4, 1); v = 0;
for g1 = 1:2
    for g2 = 1:2
        v = v + 1;
        xd = G1{g1}*Fh + xb;
        zd = zb + G2{g2}*Eh;
        switch coords
            case 'xy'
                Dphi = [2*P.phi(:,1)*x(1) + P.phi(:,2)*x(2), P.phi(:,2)*x(1) + 2*P.phi(:,3)*x(2)];
                wd = zd + Dphi*xd;
            case 'xz'
                wd = zd;
            case 'xz4'
                dm4 = [4*x(1)^3, 0; 3*x(1)^2*x(2), x(1)^3; 2*x(1)*x(2)^2, 2*x(1)^2*x(2); ...
                       x(2)^3, 3*x(1)*x(2)^2; 0, 4*x(2)^3];
                wd = zd - P.eta4*dm4*xd;
        end
        f{v} = [xd; wd];
    end
end
end

function dg = field_degree(f)
dg = 0;
for v = 1:numel(f)
    for j = 1:numel(f{v}), dg = max(dg, degree(f{v}(j))); end
end
end

function out = cs_free(P, coords, D, r, ops, enc)
[vars, f, x] = cs_system(P, coords, enc);
mons = monolist(vars, D);
dg = zeros(length(mons), 1);
for i = 1:length(mons), dg(i) = degree(mons(i)); end
mons = mons(dg >= 2);
c = sdpvar(length(mons), 1);
W = c'*mons;
ds = 2*floor((D - 1 + field_degree(f) - 2)/2);
Fc = sos(W); dec = c;
for v = 1:numel(f)
    [s, cs] = polynomial(vars, ds);
    Fc = [Fc, sos(-jacobian(W, vars)*f{v} - 1e-3*(x(1)^8 + x(2)^8) - s*(r^2 - x'*x)), sos(s)];
    dec = [dec; cs];
end
t0 = tic;
[sol, vb] = solvesos(Fc, [], ops, dec);
out.time = toc(t0); out.info = sol.info; out.problem = sol.problem; out.ok = sol.problem == 0;
out.E = monexp(mons, vars); out.c = value(c);
out.Eb = cell(1, 3); out.Eb{1} = monexp(vb{1}, vars);
Tb = []; Sb = [];
for v = 1:numel(f)
    Tb = [Tb; monexp(vb{2*v}, vars)];
    Sb = [Sb; monexp(vb{2*v + 1}, vars)];
end
out.Eb{2} = unique(Tb, 'rows'); out.Eb{3} = unique(Sb, 'rows');
end

function out = cs_manual(P, coords, E, r, central, ops, enc, opts)
if nargin < 8, opts = struct(); end
if ~isfield(opts, 'reduce'), opts.reduce = true; end
if ~isfield(opts, 'composite'), opts.composite = false; end
if ~isfield(opts, 'maxit'), opts.maxit = 10; end
[vars, f, x] = cs_system(P, coords, enc);
nv = numel(f);
% Adapted coordinates make the prescribed jet face a monomial deletion.
pc = @(Eb, dmax) all(Eb(:,3:end) == 0, 2) & sum(Eb(:,1:2), 2) <= dmax;
halfD = max(sum(E{1}, 2));
for b = 2:3
    E{b} = E{b}(sum(E{b}(:,3:end), 2) <= halfD, :);
    if opts.reduce, E{b} = E{b}(~pc(E{b}, 2), :); end
end
if opts.composite
    E{1} = E{1}(all(E{1}(:,3:end) == 0, 2) | (sum(E{1}(:,3:end), 2) == 1 & sum(E{1}(:,1:2), 2) == 0), :);
end
out.deleted = 0;
for it = 1:opts.maxit
    v1 = monvec(vars, E{1}); v2 = monvec(vars, E{2}); v3 = monvec(vars, E{3});
    n = [length(v1), length(v2), length(v3)];
    R1 = sdpvar(n(1));
    R2 = cell(nv, 1); R3 = cell(nv, 1);
    for v = 1:nv, R2{v} = sdpvar(n(2)); R3{v} = sdpvar(n(3)); end
    W = v1'*R1*v1;
    if central, e = sdpvar(1); t = sdpvar(1); else, e = 1e-3; t = 0; end
    F = [];
    if opts.composite
        cen = all(E{1}(:,3:end) == 0, 2);
        F = [F, R1(cen, ~cen) == 0];
    end
    tr = trace(R1);
    for v = 1:nv
        s = v3'*R3{v}*v3;
        T = -jacobian(W, vars)*f{v} - e*(x(1)^8 + x(2)^8) - s*(r^2 - x'*x) - v2'*R2{v}*v2;
        F = [F, coefficients(T, vars) == 0];
        tr = tr + trace(R2{v}) + trace(R3{v});
    end
    % A shared Lyapunov Gram block must satisfy every vertex; trace fixes its scale.
    if central
        F = [F, e >= t, tr == 1, R1 - t*eye(n(1)) >= 0];
        for v = 1:nv, F = [F, R2{v} - t*eye(n(2)) >= 0, R3{v} - t*eye(n(3)) >= 0]; end
        obj = -t;
    else
        F = [F, R1 >= 0];
        for v = 1:nv, F = [F, R2{v} >= 0, R3{v} >= 0]; end
        obj = [];
    end
    t0 = tic; sol = optimize(F, obj, ops); out.time = toc(t0);
    out.ok = sol.problem == 0; out.info = sol.info; out.problem = sol.problem;
    out.Q1 = value(R1); out.Q2 = cellfun(@value, R2, 'uni', 0); out.Q3 = cellfun(@value, R3, 'uni', 0);
    if central
        out.t = value(t); out.e = value(e);
        if ~out.ok || out.t > 1e-7 || it == opts.maxit, break; end
        % Optional small-diagonal deletions retry centering after the prescribed face.
        d2 = max(cell2mat(cellfun(@(Q) diag(Q)', out.Q2, 'uni', 0)), [], 1);
        d3 = max(cell2mat(cellfun(@(Q) diag(Q)', out.Q3, 'uni', 0)), [], 1);
        d1 = diag(out.Q1)';
        del = {d1 < 1e-8*max(d1), d2 < 1e-8*max(d2), d3 < 1e-8*max(d3)};
        nd = 0;
        for b = 1:3, E{b} = E{b}(~del{b}, :); nd = nd + nnz(del{b}); end
        out.deleted = out.deleted + nd;
        if nd == 0, break; end
    else
        break
    end
end
out.iter = it; out.E = E; out.n = n;
Wn = v1'*out.Q1*v1;
[cw, mw] = coefficients(Wn, vars);
out.E_W = monexp(mw, vars); out.c = value(cw);
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
