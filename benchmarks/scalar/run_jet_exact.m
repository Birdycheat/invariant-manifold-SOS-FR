% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if ~exist('seeds','var'), seeds = 1:20; end
cm_storage_guard(pwd);
addpath(fileparts(which('yalmip')), '-begin');
if ~exist('cfgs', 'var'), cfgs = [1 5; 1 7; 3 5; 3 7]; end
if ~exist('outdir', 'var'), outdir = fullfile(fileparts(pwd), 'out'); end
cm_storage_guard(outdir);
if ~exist('logf', 'var'), logf = 'jet_center_log.txt'; end
jetdir = fullfile(fileparts(pwd), 'data', 'jet_inputs');
if ~exist(jetdir, 'dir'), mkdir(jetdir); end
logp = fullfile(outdir, logf);
if exist(logp, 'file'), delete(logp); end
diary(logp);
for cfg = cfgs'
    N = cfg(1); m = cfg(2); as = 0.5*(N == 1) + 1.0*(N ~= 1);
    Bt = load(fullfile(outdir, sprintf('batch_exact_N%d_m%d.mat', N, m)));
    for seed = seeds
        row = Bt.rows(strcmp(Bt.rows(:,4), 'CM') & [Bt.rows{:,3}]' == seed, :);
        yalmip('clear');
        P = gen_instance_exact(N, m, 1, seed, [], 1, as);
        fn = fullfile(jetdir, sprintf('N%d_m%d_s%02d.json', N, m, seed));
        % Reuse the same instance and radius when centering the reduced jet program.
        [ok, out] = jet_sos_center(P, 4, row{8}, (m-1)/2, fn);
        fprintf('N=%d m=%d seed=%2d r=%.4f : ok=%d t=%.2e e=%.2e tailles=%s iter=%d %.1fs %s\n', ...
                N, m, seed, out.r, ok, out.t, out.e, mat2str(out.n), out.iter, out.time, ...
                strjoin(out.deleted', ' '));
    end
end
diary off
