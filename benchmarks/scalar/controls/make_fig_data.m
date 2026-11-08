% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(pwd);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(fileparts(pwd), 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
if ~exist('outdir', 'var'), outdir = fullfile(fileparts(pwd), 'out'); end
cm_storage_guard(outdir);
out = fullfile(fileparts(pwd), 'data', 'figdata');
if ~exist(out, 'dir'), mkdir(out); end
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',300);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr, ...
                           'frlib.reduce','auto','verbose',0);

S = load(fullfile(outdir, 'example1.mat'));
P = S.P; r = S.r;
fid = fopen(fullfile(out, 'example1_meta.txt'), 'w');
fprintf(fid, 'r %.17g\nrmax %.17g\nmu %.17g\nmumax %.17g\n', r, S.lo, S.mu, S.mumax);
fprintf(fid, 'h %s\n', char(P.HS));
fclose(fid);
yalmip('clear');
for D = [2 4 6]
    [ok, ~, info] = gen_sos_cert(P, D, r, ops0);
    wcert(out, sprintf('SOS_D%d', D), ok, info.E, info.cval, 0);
end
[ok, ~, info] = gen_sos_cert(P, 4, r, mkfr('d'));
wcert(out, 'FRd_D4', ok, info.E, info.cval, 0);
[ok, ~, info] = jp_cert(P, 6, r, ops0);
wcert(out, 'JP', ok, info.E, info.cval, 0);
for D = [4 6]
    [ok, info] = gen_sos_cert_z(P, D, r, mkfr('d'), 3);
    wcert(out, sprintf('jet_D%d', D), ok, info.E, info.cval, 3);
end

% Recompute numerical kernel-residual ratios; these diagnostics are not exact proofs.
fid = fopen(fullfile(out, 'gram_kernel.csv'), 'w');
fprintf(fid, 'N,m,seed,i,rho\n');
for cfg = [1 5; 1 7; 3 5; 3 7]'
    N = cfg(1); m = cfg(2); as = 0.5*(N == 1) + 1.0*(N ~= 1);
    B = load(fullfile(outdir, sprintf('batch_exact_N%d_m%d.mat', N, m)));
    for seed = 1:20
        row = B.rows(strcmp(B.rows(:,4), 'CM') & [B.rows{:,3}]' == seed, :);
        r = row{8};
        yalmip('clear');
        P = gen_instance_exact(N, m, 1, seed, [], 1, as);
        [ok, ~, info] = gen_sos_cert(P, 4, r, ops0);
        if ~ok, continue; end
        x = info.x; y = info.y; v = info.v{2}; Q = info.Q{2};
        Ev = zeros(length(v), 1+N);
        for k = 1:length(v)
            try
                Ev(k,1) = degree(v(k), x);
                for j = 1:N, Ev(k,1+j) = degree(v(k), y(j)); end
            catch
            end
        end
        C = jetvec(Ev, P.H, (m+1)/2);
        for i = 0:(m+1)/2
            ci = C(:, i+1);
            if norm(ci) > 0, fprintf(fid, '%d,%d,%d,%d,%.6e\n', N, m, seed, i, norm(Q*ci)/(norm(Q)*norm(ci))); end
        end
        % Compare the forced jet directions with a fixed-seed random reference.
        R = randn(RandStream('mt19937ar','Seed',seed), size(Q,1), 200);
        rho = vecnorm(Q*R)./(norm(Q)*vecnorm(R));
        fprintf(fid, '%d,%d,%d,-1,%.6e\n', N, m, seed, median(rho));   % i = -1 : aleatoire (mediane)
        fprintf('N=%d m=%d seed=%2d : ok\n', N, m, seed);
    end
end
fclose(fid);

function wcert(out, name, ok, E, c, p)
fid = fopen(fullfile(out, ['cert_' name '.csv']), 'w');
fprintf(fid, '# ok=%d p=%d\n', ok, p);
for k = 1:size(E,1), fprintf(fid, '%d,%d,%.17g\n', E(k,1), E(k,2), c(k)); end
fclose(fid);
fprintf('%s : solveur ok = %d, %d monomes\n', name, ok, size(E,1));
end

function C = jetvec(Ev, Hs, imax)
[n, nv] = size(Ev); N = nv - 1;
L = imax + 1;
hser = zeros(N, L);
for k = 1:min(imax, size(Hs,2)), hser(:, k+1) = Hs(:, k); end
C = zeros(n, L);
for i = 1:n
    s = zeros(1, L); s(1) = 1;
    ex = Ev(i,1);
    if ex > imax, continue; end
    s = [zeros(1, ex) s(1:L-ex)];
    for j = 1:N
        for e = 1:Ev(i, 1+j)
            t = conv(s, hser(j,:)); s = t(1:L);
        end
    end
    C(i,:) = s;
end
end
