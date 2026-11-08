% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(pwd);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(fileparts(pwd), 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
if ~exist('Nsel', 'var'), Nsel = 1; end
if ~exist('msel', 'var'), msel = 5; end
if ~exist('seeds', 'var'), seeds = 1:20; end
if ~exist('outdir', 'var'), outdir = fullfile(fileparts(pwd), 'out'); end
cm_storage_guard(outdir);
if ~exist(outdir, 'dir'), mkdir(outdir); end
N = Nsel; m = msel;
as = 0.5*(N == 1) + 1.0*(N ~= 1);
tag = sprintf('N%d_m%d', N, m);
jsondir = fullfile(fileparts(pwd), 'data', 'benchmark');
if ~exist(jsondir, 'dir'), mkdir(jsondir); end
logf = fullfile(outdir, sprintf('batch_exact_%s_log.txt', tag));
if exist(logf, 'file'), delete(logf); end
diary(logf);
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',300, ...
                   'mosek.MSK_IPAR_NUM_THREADS',4);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr, ...
                           'frlib.reduce','auto','verbose',0);
rows = {};
for seed = seeds
    yalmip('clear');
    P = gen_instance_exact(N, m, 1, seed, [], 1, as);
    % Bisect the anchored-family radius, then use half of it for every comparison.
    lo = 0; hi = 1.5;
    for it = 1:8
        rr = (lo+hi)/2;
        if gen_cm_cert(P, m+1, rr, ops0), lo = rr; else, hi = rr; end
    end
    if lo == 0
        fprintf('N=%d m=%d seed=%2d : aucun certificat CM, instance ignoree\n', N, m, seed);
        continue
    end
    r = 0.5*lo;
    line = sprintf('N=%d m=%d seed=%2d h2=%s r=%.3f', N, m, seed, mat2str(P.H(:,2)', 2), r);
    % Retain exact system data and the binary64 candidates for independent checks.
    J = instance_json(P, seed, r);
    [ok, mu, info] = gen_cm_cert(P, m+1, r, ops0);
    st = 'echec solveur'; vi = struct('mu', NaN, 'r', NaN);
    if ok
        [ver, vi] = cm_verify_exact(P, mu, m+1, r);
        if ver, st = 'verifie'; else, st = 'non verifie'; end
    end
    [rows, line] = addrow(rows, line, N, m, seed, 'CM', ok, st, 0, r, info.time);
    J.cm = struct('ok', ok, 'status', st, 'p', m+1, 'mu', char(sym(vi.mu)), 'r', char(sym(vi.r)));
    % Compare free and partial-FR searches in the original candidate coordinates.
    for D = [2 4 6]
        [ok, ~, info] = gen_sos_cert(P, D, r, ops0);
        [st, rb, cand] = classify_candidate(ok, info, P, r, 0);
        [rows, line] = addrow(rows, line, N, m, seed, sprintf('SOS D=%d', D), ok, st, rb, r, info.time);
        J.candidates{end+1} = candrec(sprintf('SOS D=%d', D), ok, st, rb, cand, 0);
    end
    [ok, ~, info] = jp_cert(P, 6, r, ops0);
    [st, rb, cand] = classify_candidate(ok, info, P, r, 0);
    [rows, line] = addrow(rows, line, N, m, seed, 'Jones-Peet', ok, st, rb, r, info.time);
    J.candidates{end+1} = candrec('Jones-Peet', ok, st, rb, cand, 0);
    for appr = {'d', 'dd'}
        [ok, ~, info] = gen_sos_cert(P, 4, r, mkfr(appr{1}));
        [st, rb, cand] = classify_candidate(ok, info, P, r, 0);
        [rows, line] = addrow(rows, line, N, m, seed, ['frlib ' appr{1}], ok, st, rb, r, info.time);
        J.candidates{end+1} = candrec(['frlib ' appr{1}], ok, st, rb, cand, 0);
    end
    % Hold the adapted candidate space fixed in the with/without-FR ablation.
    p = (m-1)/2;
    for v = {{'jet sans RF', ops0}, {'RF-jet', mkfr('d')}}
        [ok, info] = gen_sos_cert_z(P, 4, r, v{1}{2}, p);
        [st, rb, cand] = classify_candidate(ok, info, P, r, p);
        [rows, line] = addrow(rows, line, N, m, seed, v{1}{1}, ok, st, rb, r, info.time);
        J.candidates{end+1} = candrec(v{1}{1}, ok, st, rb, cand, p);
    end
    fprintf('%s\n', line);
    save(fullfile(outdir, sprintf('batch_exact_%s.mat', tag)), 'rows');
    fid = fopen(fullfile(jsondir, sprintf('N%d_m%d_s%02d.json', N, m, seed)), 'w');
    fprintf(fid, '%s', jsonencode(J)); fclose(fid);
end
diary off

function [rows, line] = addrow(rows, line, N, m, seed, meth, ok, st, rb, r, t)
rows(end+1,:) = {N, m, seed, meth, ok, st, rb, r, t};
ab = containers.Map({'verifie','non verifie','invalide','aucune violation','echec solveur'}, ...
                    {'V','NV','INV','ok?','X'});
line = sprintf('%s | %s %s', line, meth, ab(st));
if strcmp(st, 'invalide'), line = sprintf('%s(%.0e)', line, rb); end
end

function J = instance_json(P, seed, r)
[Lam, B, a, v, w] = exact_data(P);
q = @(A) struct('dims', size(A), 'vals', {arrayfun(@(t) char(t), A(:)', 'uni', 0)});
J = struct('N', P.N, 'm', P.m, 'seed', seed, 'r', sprintf('%.17g', r), ...
           'Lam', q(Lam), 'B', q(B), 'a', q(a), 'v', q(v), 'w', q(w), 'g', q(P.gS), ...
           'H_double', dflat(P.H), 'H_exact', q(P.HS(:, 1:min(size(P.HS,2), 12))));
J.candidates = {};
end

function s = dflat(A)
s = struct('dims', size(A), 'vals', {arrayfun(@(t) sprintf('%.17g', t), A(:)', 'uni', 0)});
end

function c = candrec(name, ok, st, rb, cand, p)
c = struct('method', name, 'solver_ok', ok, 'status', st, 'r_bad', rb, 'p', p, ...
           'E', dflat(cand.E), 'c', dflat(cand.c), ...
           'witness_x', sprintf('%.17g', cand.x), 'witness_Y', dflat(cand.Y));
end
