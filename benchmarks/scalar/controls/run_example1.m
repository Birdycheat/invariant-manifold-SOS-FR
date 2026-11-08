% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(pwd);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(fileparts(pwd), 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
if ~exist('outdir', 'var'), outdir = fullfile(fileparts(pwd), 'out'); end
cm_storage_guard(outdir);
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fullfile(outdir, 'example1_log.txt');
if exist(logf, 'file'), delete(logf); end
diary(logf);
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',300);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr, ...
                           'frlib.reduce','auto','verbose',0);
K = 24; m = 7;
P = struct('N',1,'m',m,'d',1,'seed',0,'lam',1,'Lam',-1,'B',0,'a',1,'v',2,'w',1, ...
           'g',[0 0 0 -2 -1 0 0 -1],'H',zeros(1,K));
P.LamS = sym(-1); P.BS = sym(0); P.aS = sym(1); P.vS = sym(2); P.wS = sym(1);
[P.gS, P.HS] = exact_jet(P, K, false);
P.H = double(P.HS);
fprintf('h_2..h_10 = %s\n', char(P.HS(2:10)));

% Search the anchored radius, then certify the rounded parameters by exact root counts.
pcm = m - 1;
lo = 0; hi = 2;
for it = 1:12
    rr = (lo+hi)/2;
    if gen_cm_cert(P, pcm, rr, ops0), lo = rr; else, hi = rr; end
end
[~, mumax] = gen_cm_cert(P, pcm, lo, ops0);
[vermax, vimax] = cm_verify_exact(P, mumax, pcm, lo);
fprintf('CM : r_max = %.4f (mu = %.4f), verification exacte (Sturm) : %d en %.1f s\n', lo, mumax, vermax, vimax.time);
xs = linspace(-2, 2, 40001); best = 0;
for mu = linspace(0.126, 2, 2000)
    bad = 8*mu*(1+4*xs.^2+3*xs.^3) - (1+8*mu*xs.^2+6*mu*xs.^3).^2 <= 0;
    best = max(best, min(abs(xs(bad))));
end
fprintf('forme close (sans marge) : sup_mu r(mu) ~ %.4f\n', best);
r = 0.5*lo;
fprintf('bande commune : r = %.4f\n\n', r);

% Use the same prescribed half-radius across all numerical formulations.
rows = {};
[ok, mu, info] = gen_cm_cert(P, pcm, r, ops0);
[ver, vi] = cm_verify_exact(P, mu, pcm, r);
rows(end+1,:) = {'CM (Thm 4)', ok, tern(ver, 'verifie (Sturm)', 'NON verifie'), 0, info.time};
for D = [2 4 6]
    [ok, ~, info] = gen_sos_cert(P, D, r, ops0);
    [st, rb] = classify_candidate(ok, info, P, r, 0);
    rows(end+1,:) = {sprintf('SOS libre D=%d', D), ok, st, rb, info.time}; %#ok<SAGROW>
end
for appr = {'d', 'dd'}
    for D = [4 6]
        [ok, ~, info] = gen_sos_cert(P, D, r, mkfr(appr{1}));
        [st, rb] = classify_candidate(ok, info, P, r, 0);
        rows(end+1,:) = {sprintf('frlib %s, D=%d', appr{1}, D), ok, st, rb, info.time}; %#ok<SAGROW>
    end
end
[ok, ~, info] = jp_cert(P, 6, r, ops0);
[st, rb] = classify_candidate(ok, info, P, r, 0);
rows(end+1,:) = {'Jones-Peet', ok, st, rb, info.time};
p = (m-1)/2;
for D = [4 6]
    [ok, info] = gen_sos_cert_z(P, D, r, mkfr('d'), p);
    [st, rb] = classify_candidate(ok, info, P, r, p);
    rows(end+1,:) = {sprintf('RF-jet (z, p=3), D=%d', D), ok, st, rb, info.time}; %#ok<SAGROW>
end
fprintf('%-24s %-8s %-20s %-10s %s\n', 'methode', 'solveur', 'statut exact', 'r_bad', 'temps (s)');
for i = 1:size(rows,1)
    fprintf('%-24s %-8d %-20s %-10.1e %.2f\n', rows{i,:});
end
jetdir = fullfile(fileparts(pwd), 'data', 'jet_inputs');
for D = [4 6]
    yalmip('clear');
    [ok, out] = jet_sos_center(P, D, r, p, fullfile(jetdir, sprintf('ex1_D%d.json', D)));
    fprintf('SOS par jet, point central, D=%d : ok=%d t=%.2e tailles=%s\n', D, ok, out.t, mat2str(out.n));
end
save(fullfile(outdir, 'example1.mat'), 'P', 'lo', 'mumax', 'r', 'mu', 'rows');
diary off

function o = tern(c, a, b)
if c, o = a; else, o = b; end
end
