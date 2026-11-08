% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(outdir);
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(here, '..', 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
instdir = fullfile(here, '..', 'data', 'triple', 'instances');
if ~exist('files', 'var')
    L = dir(fullfile(instdir, 'tr_ny*_s*.json')); files = {L.name};
end
canddir = fullfile(outdir, 'candidates');
if ~exist(canddir, 'dir'), mkdir(canddir); end
diary(fullfile(outdir, 'triple_batch_log.txt'));
fprintf('Output: %s\n', outdir);
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',1200, ...
                   'mosek.MSK_IPAR_NUM_THREADS',4);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr, ...
                           'frlib.reduce','auto','verbose',0);
D = 4;
for fi = 1:numel(files)
    yalmip('clear');
    S = jsondecode(fileread(fullfile(instdir, files{fi})));
    P = numdata(S);
    r = P.r;
    if exist(fullfile(canddir, files{fi}), 'file'), continue; end
    J = struct('instance', files{fi}, 'r', S.r, 'D', D, 'class', S.leading.class, 'candidates', {{}});
    line = sprintf('%-18s %-9s', files{fi}, S.leading.class);
    [ok, info] = md_free(P, 'xy', D, r, ops0);
    J.candidates{end+1} = cand('SOS', 'xy', ok, info);
    Exy = bases_xy(P, info);
    [ok2, info2] = md_free(P, 'xy', D, r, mkfr('d'));
    J.candidates{end+1} = cand('frlib d', 'xy', ok2, info2);
    if P.ny == 1
        [ok2, info2] = md_manual(P, 'xy', Exy, r, false, ops0, 'jet');
        J.candidates{end+1} = cand('Thm2', 'xy', ok2, info2);
    end
    [ok, info] = md_free(P, 'xz', D, r, ops0);
    J.candidates{end+1} = cand('jet', 'xz', ok, info);
    Exz = bases_xz(P);
    [ok2, info2] = md_manual(P, 'xz', Exz, r, false, ops0, 'jet');
    J.candidates{end+1} = cand('Thm3', 'xz', ok2, info2);
    J.t = struct(); J.t_solver = struct();
    runs = {{'xz', Exz, 'none', 'xz_none'}};
    if P.ny == 1
        runs = [runs, {{'xy', Exy, 'none', 'xy_none'}, {'xy', Exy, 'jet', 'xy_thm2'}}];
    end
    for v = runs
        [~, info3] = md_manual(P, v{1}{1}, v{1}{2}, r, true, ops0, v{1}{3}, 1);
        J.t.(v{1}{4}) = info3.t; J.t_solver.(v{1}{4}) = info3.info;
    end
    [ok2, info2] = md_manual(P, 'xz', Exz, r, true, ops0, 'jet');
    J.t.xz_thm3 = info2.t; J.t.xz_thm3_first = info2.t_first;
    % Export the bases actually used together with their matching Gram matrices.
    J.center = struct('ok', ok2, 't', info2.t, 'e', info2.e, 'iter', info2.iter, ...
        'deleted', {info2.deleted_monomials}, ...
        'E1', info2.E{1}, 'E2', info2.E{2}, 'E3', info2.E{3}, ...
        'Q1', flat(info2.Q{1}), 'Q2', flat(info2.Q{2}), 'Q3', flat(info2.Q{3}), 'time', info2.time);
    for c = J.candidates, line = [line sprintf(' | %s:%d %.0fs', c{1}.name, c{1}.ok, c{1}.time)]; end %#ok<AGROW>
    if P.ny == 1, line = [line sprintf(' | xy %.1e / Thm2 %.1e', J.t.xy_none, J.t.xy_thm2)]; end
    fprintf('%s | xz %.1e / Thm3 %.1e (first %.1e, %d it.) sizes %s\n', line, ...
            J.t.xz_none, J.t.xz_thm3, J.t.xz_thm3_first, info2.iter, mat2str(info2.n));
    fid = fopen(fullfile(canddir, files{fi}), 'w'); fprintf(fid, '%s', jsonencode(J)); fclose(fid);
end
diary off

function P = numdata(S)
N = S.num; P.ny = S.ny; P.r = N.r;
P.a = reshape(N.a, 3, S.ny); P.Lam = reshape(N.Lam, S.ny, S.ny);
P.B = {reshape(N.B1, S.ny, S.ny), reshape(N.B2, S.ny, S.ny), reshape(N.B3, S.ny, S.ny)};
P.phi = reshape(N.phi, S.ny, 6); P.F = reshape(N.Fhat, 3, 21); P.E = reshape(N.Ehat, S.ny, 15);
P.L5 = N.layer5; P.L4 = N.layer4; P.L2 = N.layer2;
end

function m = monos(x, L)
m = [];
for i = 1:size(L, 1)
    t = 1;
    for j = 1:3, if L(i,j) > 0, t = t*x(j)^L(i,j); end, end
    m = [m; t];
end
end

function [vars, f, x] = md_system(P, coords)
x = sdpvar(3, 1); w = sdpvar(P.ny, 1);
phi = P.phi*monos(x, P.L2);
if strcmp(coords, 'xy'), z = w - phi; else, z = w; end
m5 = monos(x, P.L5); m4 = monos(x, P.L4);
xd = P.F*m5 + x.*(P.a*z);
zd = (P.Lam + x(1)*P.B{1} + x(2)*P.B{2} + x(3)*P.B{3})*z + P.E*m4;
% Pull back the transverse equation when using the original coordinates.
if strcmp(coords, 'xy')
    wd = zd + jacobian(phi, x)*xd;
else
    wd = zd;
end
vars = [x; w]; f = [xd; wd];
end

function E = low_exponents(nv, dmax)
E = zeros(1, nv);
for d = 1:dmax
    E = [E; compositions(nv, d)];
end
end

function C = compositions(n, d)
if n == 1, C = d; return; end
C = zeros(0, n);
for k = d:-1:0
    R = compositions(n - 1, d - k);
    C = [C; [k*ones(size(R,1),1), R]];
end
end

function E = bases_xz(P)
nv = 3 + P.ny;
E1 = low_exponents(nv, 2); E1 = E1(sum(E1, 2) >= 1, :);
E2 = low_exponents(nv, 4); E2 = E2(sum(E2(:,4:end), 2) <= 2, :);
E3 = low_exponents(nv, 3); E3 = E3(sum(E3(:,4:end), 2) <= 2, :);
E = {E1, E2, E3};
end

function E = bases_xy(P, info)
nv = 3 + P.ny;
if isfield(info, 'Eb') && numel(info.Eb) == 3
    E = info.Eb;
else
    E1 = low_exponents(nv, 2); E1 = E1(sum(E1, 2) >= 1, :);
    E = {E1, low_exponents(nv, 5), low_exponents(nv, 4)};
end
for b = 2:3, E{b} = E{b}(sum(E{b}(:,4:end), 2) <= 2, :); end
end

function [ok, info] = md_free(P, coords, D, r, ops)
[vars, f, x] = md_system(P, coords);
mons = monolist(vars, D);
dg = zeros(length(mons), 1);
for i = 1:length(mons), dg(i) = degree(mons(i)); end
mons = mons(dg >= 2);
c = sdpvar(length(mons), 1);
W = c'*mons;
dW = jacobian(W, vars)*f;
ds = 2*floor((degree(dW) - 2)/2);
[s, cs] = polynomial(vars, ds);
Fc = [sos(W), sos(-dW - 1e-3*(x(1)^8 + x(2)^8 + x(3)^8) - s*(r^2 - x'*x)), sos(s)];
t0 = tic;
[sol, vb] = solvesos(Fc, [], ops, [c; cs]);
info.time = toc(t0); info.info = sol.info; info.problem = sol.problem;
ok = sol.problem == 0;
info.E = monexp(mons, vars); info.c = value(c);
if iscell(vb), info.Eb = cellfun(@(q) monexp(q, vars), vb, 'uni', 0); end
end

function [ok, out] = md_manual(P, coords, E, r, central, ops, reduce, maxit)
if nargin < 8, maxit = 10; end
[vars, f, x] = md_system(P, coords);
pc = @(Eb, dmax) all(Eb(:,4:end) == 0, 2) & sum(Eb(:,1:3), 2) <= dmax;
% In adapted coordinates the prescribed face removes low central monomials.
if strcmp(reduce, 'jet') && strcmp(coords, 'xz')
    for b = 2:3, E{b} = E{b}(~pc(E{b}, 2), :); end
end
out.deleted = 0; out.deleted_monomials = {}; out.t_first = NaN;
for it = 1:maxit
    Eused = E;
    V = cellfun(@(Eb) eye(size(Eb, 1)), E, 'uni', 0);
    if strcmp(reduce, 'jet') && strcmp(coords, 'xy')
        for b = 2:3
            C = jetvectors(E{b}, P);
            if any(C(:)), V{b} = null(C'); end
        end
    end
    v = cellfun(@(Eb) monvec(vars, Eb), E, 'uni', 0);
    w = cellfun(@(Vb, vb) Vb'*vb, V, v, 'uni', 0);
    n = cellfun(@(Vb) size(Vb, 2), V);
    R = arrayfun(@(nb) sdpvar(nb), n, 'uni', 0);
    W = w{1}'*R{1}*w{1}; s = w{3}'*R{3}*w{3};
    if central, e = sdpvar(1); t = sdpvar(1); else, e = 1e-3; end
    T = -jacobian(W, vars)*f - e*(x(1)^8 + x(2)^8 + x(3)^8) - s*(r^2 - x'*x) - w{2}'*R{2}*w{2};
    F = coefficients(T, vars) == 0;
    if central
        % Trace normalization fixes scale before maximizing the common Gram margin.
        F = [F, e >= t, trace(R{1}) + trace(R{2}) + trace(R{3}) == 1];
        for b = 1:3, F = [F, R{b} - t*eye(n(b)) >= 0]; end
        obj = -t;
    else
        for b = 1:3, F = [F, R{b} >= 0]; end
        obj = [];
    end
    t0 = tic; sol = optimize(F, obj, ops); out.time = toc(t0);
    ok = sol.problem == 0; out.info = sol.info; out.problem = sol.problem;
    out.Q = cellfun(@(Vb, Rb) Vb*value(Rb)*Vb', V, R, 'uni', 0);
    if central
        out.t = value(t); out.e = value(e);
        if it == 1, out.t_first = out.t; end
        if ~ok || out.t > 1e-7, break; end
        nd = 0;
        for b = 1:3
            dgq = diag(out.Q{b}); del = dgq < 1e-8*max(dgq);
            out.deleted_monomials = [out.deleted_monomials; ...
                arrayfun(@(i) sprintf('Q%d:%s', b, mat2str(E{b}(i,:))), find(del)', 'uni', 0)'];
            E{b} = E{b}(~del, :); nd = nd + nnz(del);
        end
        out.deleted = out.deleted + nd;
        if nd == 0, break; end
    else
        break
    end
end
out.iter = it; out.E = Eused; out.n = n;
Wn = monvec(vars, Eused{1})'*out.Q{1}*monvec(vars, Eused{1});
[cw, mw] = coefficients(Wn, vars);
out.E_W = monexp(mw, vars); out.c = value(cw);
end

function C = jetvectors(Eb, P)
alphas = [low_exponents(3, 2)];
C = zeros(size(Eb, 1), size(alphas, 1));
L2 = P.L2;
for i = 1:size(Eb, 1)
    ex = Eb(i, 1:3); et = Eb(i, 4:end);
    if all(et == 0)
        [tf, k] = ismember(ex, alphas, 'rows');
        if tf, C(i, k) = 1; end
    elseif sum(et) == 1 && all(ex == 0)
        for j = 1:size(L2, 1)
            [~, k] = ismember(L2(j,:), alphas, 'rows');
            C(i, k) = P.phi(et == 1, j);
        end
    end
end
C = C(:, any(C, 1));
end

function c = cand(name, coords, ok, info)
c = struct('name', name, 'coords', coords, 'ok', ok, 'problem', info.problem, 'info', info.info, 'time', info.time);
if isfield(info, 'E_W'), E = info.E_W; else, E = info.E; end
c.E = E; c.c = arrayfun(@(q) sprintf('%.17g', q), info.c(:)', 'uni', 0);
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

function s = flat(A)
s = struct('dims', size(A), 'vals', {arrayfun(@(q) sprintf('%.17g', q), A(:)', 'uni', 0)});
end
