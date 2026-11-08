% Native MATLAB distribution: historical numerical algorithm, MIT project license.
cm_storage_guard(outdir);
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fileparts(which('yalmip')), '-begin');
frdir = fullfile(here, '..', 'third_party', 'frlib');
addpath(frdir, fullfile(frdir, 'util'));
if ~exist('instdir', 'var'), instdir = fullfile(here, '..', 'data', 'coupled', 'instances'); end
if ~exist('enc', 'var'), enc = 1; end
if ~exist('files', 'var')
    L = dir(fullfile(instdir, 'cp_ny*_s*.json')); files = {L.name};
end
canddir = fullfile(outdir, 'candidates'); centdir = fullfile(outdir, 'center');
if ~exist(canddir, 'dir'), mkdir(canddir); end
if ~exist(centdir, 'dir'), mkdir(centdir); end
diary(fullfile(outdir, 'coupled_batch_log.txt'));
fprintf('Output: %s ; enc = %d\n', outdir, enc);
ops0 = sdpsettings('solver','mosek','verbose',0,'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',300, ...
                   'mosek.MSK_IPAR_NUM_THREADS',4);
mkfr = @(appr) sdpsettings('solver','frlib','frlib.solver','mosek','frlib.approximation',appr, ...
                           'frlib.reduce','auto','verbose',0);
% Keep numerical feasibility comparisons separate from central-point reconstruction.
D = 4;
for fi = 1:numel(files)
    S = jsondecode(fileread(fullfile(instdir, files{fi})));
    P = coupled_sos('data', S);
    r = P.r;
    J = struct('instance', files{fi}, 'r', S.r, 'D', D, 'enc', enc, 'class', S.leading.class, 'candidates', {{}});
    line = sprintf('%-18s %-9s', files{fi}, S.leading.class);
    yalmip('clear');
    o = coupled_sos('free', P, 'xy', D, r, ops0, enc);
    J.candidates{end+1} = cand('SOS', 'xy', o);
    Exy = o.Eb;
    for appr = {'d', 'dd'}
        yalmip('clear');
        o = coupled_sos('free', P, 'xy', D, r, mkfr(appr{1}), enc);
        J.candidates{end+1} = cand(['frlib ' appr{1}], 'xy', o);
    end
    yalmip('clear');
    o = coupled_sos('free', P, 'xz', D, r, ops0, enc);
    J.candidates{end+1} = cand('jet', 'xz', o);
    Exz = o.Eb;
    yalmip('clear');
    o = coupled_sos('manual', P, 'xz', Exz, r, false, ops0, enc);
    J.candidates{end+1} = cand('Thm3', 'xz', o);
    yalmip('clear');
    % Compare reduced and unreduced centering within the same adapted basis family.
    none = coupled_sos('manual', P, 'xz', Exz, r, true, ops0, enc, struct('reduce', false, 'maxit', 1));
    yalmip('clear');
    none_xy = coupled_sos('manual', P, 'xy', Exy, r, true, ops0, enc, struct('reduce', false, 'maxit', 1));
    J.t = struct('xz_none', none.t, 'xy_none', none_xy.t);
    yalmip('clear');
    c = coupled_sos('manual', P, 'xz', Exz, r, true, ops0, enc);
    J.t.xz_thm3 = c.t; J.t.iter = c.iter; J.t.deleted = c.deleted; J.t.center_time = c.time;
    % Retain the central point and its ordered Gram bases for exact reconstruction.
    write_center(fullfile(centdir, files{fi}), files{fi}, 'xz', enc, S.r, c, false);
    for k = 1:numel(J.candidates)
        q = J.candidates{k}; line = [line sprintf(' | %s:%d %.1fs', q.name, q.ok, q.time)]; %#ok<SAGROW>
    end
    fprintf('%s | t: xy %.1e, xz %.1e / Thm3 %.1e (%d it.)\n', line, J.t.xy_none, J.t.xz_none, J.t.xz_thm3, c.iter);
    fid = fopen(fullfile(canddir, files{fi}), 'w'); fprintf(fid, '%s', jsonencode(J)); fclose(fid);
end
diary off

function c = cand(name, coords, info)
c = struct('name', name, 'coords', coords, 'ok', info.ok, 'problem', info.problem, ...
           'info', info.info, 'time', info.time);
if isfield(info, 'E_W'), E = info.E_W; else, E = info.E; end
c.E = E; c.c = arrayfun(@(q) sprintf('%.17g', q), info.c(:)', 'uni', 0);
end
