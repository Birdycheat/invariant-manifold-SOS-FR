% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(outdir);
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(here, '..', 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
instdir = fullfile(here, '..', 'data', 'beg', 'instances');
if ~exist('files', 'var')
    L = dir(fullfile(instdir, 'a_*.json')); files = {L.name};
    files = ['a_2_1.json', setdiff(files, {'a_2_1.json'}, 'stable')];
end
logf = fullfile(outdir, 'beg_batch_log.txt');
fprintf('Output: %s\n', outdir);
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',600,'mosek.MSK_IPAR_NUM_THREADS',4);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr,'frlib.reduce','auto','verbose',0);
for fi = 1:numel(files)
    FD = jsondecode(fileread(fullfile(instdir, files{fi})));
    [~, stem] = fileparts(files{fi});
    cdir = fullfile(outdir, stem);
    if ~exist(cdir, 'dir'), mkdir(cdir); end
    note(logf, sprintf('=== %s (%s)', stem, FD.class));
    r = 1/8;
    for D = [2 4 6], run_free(FD, 'xy', D, r, [], ops0, sprintf('A_SOS_D%d', D), cdir, logf); end
    run_jp(FD, 6, r, ops0, 'A_Jones_Peet', cdir, logf);
    run_free(FD, 'xy', 4, r, [], mkfr('d'), 'A_frlib_d', cdir, logf);
    run_free(FD, 'xy', 4, r, [], mkfr('dd'), 'A_frlib_dd', cdir, logf);
    run_free(FD, 'xz', 4, r, [], ops0, 'A_jet_noFR', cdir, logf);
    run_free(FD, 'xz', 4, r, [], mkfr('d'), 'A_jet_FRd', cdir, logf);
    B0 = bases_of(cdir, 'A_jet_noFR', outdir);
    run_center(FD, 'xz', B0, 2, r, [], ops0, 'A_center_q2', cdir, logf);
    run_center(FD, 'xz', B0, -1, r, [], ops0, 'A_center_none', cdir, logf);
    % Replace the surrogate by four polynomial vertices enclosing the true field.
    radii = FD.radii;
    r16 = radii(1);
    run_free(FD, 'vertices', 4, r16.rf, r16.Mf, ops0, 'B_jet_noFR', cdir, logf);
    run_free(FD, 'vertices', 4, r16.rf, r16.Mf, mkfr('d'), 'B_jet_FRd', cdir, logf);
    B1 = bases_of(cdir, 'B_jet_noFR', outdir);
    run_center(FD, 'vertices', B1, -1, r16.rf, r16.Mf, ops0, 'B_center_none_r16', cdir, logf);
    for ir = 1:numel(radii)
        run_center(FD, 'vertices', B1, 2, radii(ir).rf, radii(ir).Mf, ops0, sprintf('B_center_q2_r%d', round(1/radii(ir).rf)), cdir, logf);
    end
end

function P = poly_of(rows, comp, x, y)
P = 0*x;
for i = 1:numel(rows)
    e = rows{i};
    if e{1} == comp, P = P + e{5}*x^e{2}*y^e{3}; end
end
end

function [vars, fs, x] = system_of(FD, coords, M)
x = sdpvar(1); w = sdpvar(1); vars = [x; w];
switch coords
    case 'xy', fs = {[poly_of(FD.Sxy, 1, x, w); poly_of(FD.Sxy, 2, x, w)]};
    case 'xz', fs = {[poly_of(FD.Sxz, 1, x, w); poly_of(FD.Sxz, 2, x, w)]};
    case 'vertices'
        Px = poly_of(FD.T5xz, 1, x, w); Pz = poly_of(FD.T5xz, 2, x, w);
        rho3 = (x^2 + w^2)^3;
        fs = {};
        for s1 = [-1 1]
            for s2 = [-1 1]
                fs{end+1} = [Px + s1*M*rho3; Pz + s2*M*rho3];
            end
        end
end
end

function run_free(FD, coords, D, r, M, ops, name, cdir, logf)
yalmip('clear');
[vars, fs, x] = system_of(FD, coords, M);
mons = monolist(vars, D);
dg = zeros(length(mons), 1);
for i = 1:length(mons), dg(i) = degree(mons(i)); end
mons = mons(dg >= 2);
c = sdpvar(length(mons), 1);
W = c'*mons;
k = min(8, D + 4);
gam = r^2 - vars'*vars;
Fc = sos(W); prm = c;
for v = 1:numel(fs)
    dW = jacobian(W, vars)*fs{v};
    ds = 2*floor((degree(dW) - 2)/2);
    [s, cs] = polynomial(vars, ds);
    Fc = [Fc, sos(-dW - 1e-3*x^k - s*gam), sos(s)];
    prm = [prm; cs];
end
t0 = tic;
[sol, vb] = solvesos(Fc, [], ops, prm);
R = struct('name', name, 'coords', coords, 'D', D, 'k', k, 'r', r, 'problem', sol.problem, 'info', sol.info, 'seconds', toc(t0));
R.E = monexp(mons, vars); R.c = value(c);
if iscell(vb), R.Eb = cellfun(@(q) monexp(q, vars), vb, 'uni', 0); else, R.Eb = {}; end
write_json(fullfile(cdir, [name '.json']), R);
note(logf, sprintf('%-18s pb=%d %6.1fs | %s', name, sol.problem, R.seconds, sol.info));
end

function run_jp(FD, d2, R0, ops, name, cdir, logf)
yalmip('clear');
[vars, fs, ~] = system_of(FD, 'xy', []);
x = vars(1); y = vars(2); f = fs{1};
mx = monolist(x, d2);
cJ = sdpvar(d2, 1); J = cJ'*mx(2:end);
cH = sdpvar(1, d2); H = cH*mx(1:d2);
Pq = sdpvar(1);
V = J + y*H + y*Pq*y;
dV = jacobian(V, vars)*f;
[s, cs] = polynomial(vars, 2*floor((degree(dV) - 2)/2));
Fc = [sos(V - 1e-5*(x^2 + y^2)), sos(-dV - s*(R0^2 - x^2 - y^2)), sos(s)];
t0 = tic;
[sol, vb] = solvesos(Fc, [], ops, [cJ; cH(:); Pq; cs]);
E = zeros(0, 2); c = zeros(0, 1);
for k = 1:d2, E(end+1,:) = [k 0]; c(end+1,1) = value(cJ(k)); end
for k = 1:d2, E(end+1,:) = [k-1 1]; c(end+1,1) = value(cH(k)); end
E(end+1,:) = [0 2]; c(end+1,1) = value(Pq);
R = struct('name', name, 'coords', 'xy', 'D', d2, 'k', 0, 'r', R0, 'problem', sol.problem, 'info', sol.info, 'seconds', toc(t0), 'E', E, 'c', c);
if iscell(vb), R.Eb = cellfun(@(q) monexp(q, vars), vb, 'uni', 0); else, R.Eb = {}; end
write_json(fullfile(cdir, [name '.json']), R);
note(logf, sprintf('%-18s pb=%d %6.1fs | %s', name, sol.problem, R.seconds, sol.info));
end

function B = bases_of(cdir, name, outdir)
J = jsondecode(fileread(fullfile(cdir, [name '.json'])));
B = J.Eb;
if J.problem ~= 0 || numel(B) < 3 || isempty(B{2})
    J = jsondecode(fileread(fullfile(outdir, 'a_2_1', [name '.json'])));
    B = J.Eb;
end
end

function run_center(FD, coords, B0, q, r, M, ops, name, cdir, logf)
yalmip('clear');
[vars, fs, x] = system_of(FD, coords, M);
nv = numel(fs);
E1 = mat(B0{1}); E2 = mat(B0{2}); E3 = mat(B0{3});
% Remove only the prescribed low central monomials when q selects the jet face.
low = @(Eb) Eb(:,2) == 0 & Eb(:,1) <= q;
E2 = E2(~low(E2), :); E3 = E3(~low(E3), :);
gam = r^2 - vars'*vars;
deleted = {};
for it = 1:6
    n1 = size(E1,1); n2 = size(E2,1); n3 = size(E3,1);
    v1 = monvec(vars, E1); v2 = monvec(vars, E2); v3 = monvec(vars, E3);
    Q1 = sdpvar(n1); e = sdpvar(1); t = sdpvar(1);
    Q2 = arrayfun(@(i) sdpvar(n2), 1:nv, 'uni', 0); Q3 = arrayfun(@(i) sdpvar(n3), 1:nv, 'uni', 0);
    W = v1'*Q1*v1;
    F = [e >= t, Q1 - t*eye(n1) >= 0]; tr = trace(Q1);
    for v = 1:nv
        s = v3'*Q3{v}*v3;
        T = -jacobian(W, vars)*fs{v} - e*x^8 - s*gam - v2'*Q2{v}*v2;
        F = [F, coefficients(T, vars) == 0, Q2{v} - t*eye(n2) >= 0, Q3{v} - t*eye(n3) >= 0];
        tr = tr + trace(Q2{v}) + trace(Q3{v});
    end
    % Maximize a common margin across all vertices at fixed total Gram trace.
    F = [F, tr == 1];
    t0 = tic; sol = optimize(F, -t, ops); sec = toc(t0);
    tv = value(t);
    if it == 1, tfirst = tv; end
    Eused = {E1, E2, E3};
    if sol.problem ~= 0 || tv > 1e-7 || q < 0, break; end
    dg2 = max(cell2mat(cellfun(@(Q) diag(value(Q)), Q2, 'uni', 0)), [], 2);
    dg3 = max(cell2mat(cellfun(@(Q) diag(value(Q)), Q3, 'uni', 0)), [], 2);
    dg1 = diag(value(Q1));
    d1 = dg1 < 1e-8*max(dg1); d2 = dg2 < 1e-8*max(dg2); d3 = dg3 < 1e-8*max(dg3);
    if nnz(d1) + nnz(d2) + nnz(d3) == 0 || it == 6, break; end
    deleted = [deleted; arrayfun(@(i) sprintf('Q1:%s', mat2str(E1(i,:))), find(d1)', 'uni', 0)'; ...
        arrayfun(@(i) sprintf('Q2:%s', mat2str(E2(i,:))), find(d2)', 'uni', 0)'; ...
        arrayfun(@(i) sprintf('Q3:%s', mat2str(E3(i,:))), find(d3)', 'uni', 0)']; %#ok<AGROW>
    E1 = E1(~d1,:); E2 = E2(~d2,:); E3 = E3(~d3,:);
end
R = struct('name', name, 'coords', coords, 'q', q, 'r', r, 'M', M, 'k', 8, 'vertices', nv, 'problem', sol.problem, 'info', sol.info, 'seconds', sec, 't', tv, 't_first', tfirst, 'e', value(e), 'iter', it);
R.deleted = deleted; R.E1 = Eused{1}; R.E2 = Eused{2}; R.E3 = Eused{3};
R.Q1 = value(Q1); R.Q2 = cellfun(@value, Q2, 'uni', 0); R.Q3 = cellfun(@value, Q3, 'uni', 0);
write_json(fullfile(cdir, [name '.json']), R);
note(logf, sprintf('%-18s pb=%d t=%.3e first=%.3e e=%.3e it=%d sizes=[%d %d %d]x%d %5.1fs deleted=%s', name, sol.problem, tv, tfirst, value(e), it, size(Eused{1},1), size(Eused{2},1), size(Eused{3},1), nv, sec, strjoin(deleted', ' ')));
end

function M = mat(A)
M = zeros(numel(A), 2);
for i = 1:numel(A), M(i,:) = str2double(A{i}(:))'; end
end

function note(logf, s)
fprintf('%s\n', s);
fid = fopen(logf, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
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

% Keep full-precision numerical entries and ordered bases for exact reconstruction.
function write_json(fn, R)
fns = fieldnames(R);
for i = 1:numel(fns)
    v = R.(fns{i});
    if isnumeric(v) && ~isscalar(v) && ~isempty(v), R.(fns{i}) = num2cell_rows(v); end
    if iscell(v) && ~isempty(v) && isnumeric(v{1})
        R.(fns{i}) = cellfun(@num2cell_rows, v, 'uni', 0);
    end
end
fid = fopen(fn, 'w'); fprintf(fid, '%s', jsonencode(R, 'ConvertInfAndNaN', true)); fclose(fid);
end

function c = num2cell_rows(A)
c = cell(size(A,1), 1);
for i = 1:size(A,1), c{i} = arrayfun(@(q) sprintf('%.17g', q), A(i,:), 'uni', 0); end
end
