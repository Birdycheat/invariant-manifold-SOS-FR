% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if ~exist('seeds','var'), seeds = 1:20; end
cm_storage_guard(outdir);
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fileparts(which('yalmip')), '-begin');
jetdir = fullfile(here, '..', 'data', 'jet_inputs');
distributionRoot = fileparts(here);
if ~exist('cfgs', 'var'), cfgs = [1 5; 1 7; 3 5; 3 7]; end
if ~exist('tag', 'var'), tag = 'all'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
cd(outdir);
fprintf('Sorties : %s\n', outdir);
if ~exist('scalar_inputs', 'dir'), mkdir('scalar_inputs'); end
logf = fopen(sprintf('scalar_regions_%s.jsonl', tag), 'a');
D = 4; NBIS = 7; TMIN = 1e-8;
for cfg = cfgs'
    N = cfg(1); m = cfg(2); as = 0.5*(N == 1) + 1.0*(N ~= 1); p = (m - 1)/2;
    for seed = seeds
        name = sprintf('N%d_m%d_s%02d', N, m, seed);
        J = jsondecode(fileread(fullfile(jetdir, [name '.json'])));
        rband = str2rat(J.r);
        t0 = tic;
        yalmip('clear');
        P = gen_instance_exact(N, m, 1, seed, [], 1, as);
        lo = rband; hi = 8*rband;
        fn = fullfile('scalar_inputs', [name '.json']);
        % Numerical bisection proposes a radius from the centering margin.
        for it = 1:NBIS
            rr = (lo + hi)/2;
            yalmip('clear');
            [ok, o] = jet_sos_center(P, D, rr, p, fn);
            if ok && o.t > TMIN, lo = rr; else, hi = rr; end
        end
        proved = false; rr = lo; nround = 0; tval = NaN;
        % Back off until rational reconstruction proves the proposed enlarged strip.
        while ~proved && nround < 8 && rr > rband
            nround = nround + 1;
            yalmip('clear');
            [ok, o] = jet_sos_center(P, D, rr, p, fn);
            tval = o.t;
            if ok
                roundReport = cm_round_certificates(distributionRoot, fullfile(outdir,fn), fullfile(outdir,'certificates'), 'scalar');
                proved = roundReport.proved == 1;
            end
            if ~proved, rr = rr - 0.03*rband; end
        end
        rp = floor(1e4*rr)/1e4;
        rec = struct('instance', name, 'r_band', rband, 'r_max_cm', 2*rband, 'r_float', lo, ...
                     'r_proved', rp*proved, 'proved', proved, 't', tval, 'rounds', nround, 'seconds', toc(t0));
        fprintf(logf, '%s\n', jsonencode(rec));
        fprintf('%s r_max(CM)=%.4f r_float=%.4f r_proved=%.4f t=%.1e %.0fs\n', name, 2*rband, lo, rec.r_proved, tval, rec.seconds);
    end
end
fclose(logf);

function v = str2rat(s)
if isnumeric(s), v = s; return; end
parts = strsplit(s, '/');
v = str2double(parts{1});
if numel(parts) > 1, v = v/str2double(parts{2}); end
end
