function report = review_dd_scaling(outdir, report_path)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
here = fileparts(mfilename('fullpath'));
addpath(here);
if nargin < 1 || isempty(outdir), outdir = tempname; end
cm_storage_guard(outdir);
if ~exist(outdir, 'dir'), mkdir(outdir); end
if nargin < 2 || isempty(report_path), report_path = fullfile(outdir, 'summary.json'); end
report = struct('experiment', 'FR dd with diagonal scaling, fixed polynomial spaces', ...
    'schema_version',2, ...
    'status', 'running', 'raw_directory', outdir, 'matlab', version, ...
    'candidate_coordinates', 'original x,y', 'candidate_degree', 4, ...
    'margin_weight', '1/1000', 'normalization', 'none, as in gen_sos_cert', ...
    'scaling_rule', 's_i=1/max_d |C(i,d)| if nonzero, otherwise 1', ...
    'sdp_orientation', 'YALMIP optimize dualize=0; frlib auto therefore reduces the dual', ...
    'interpretation', 'This explicit SDP is not a replay of the archived solvesos baseline.', ...
    'timing_scope', 'optimize call only; excludes basis extraction, symbolic jet setup and witness classification', ...
    'scope', 'Two representative scalar instances; not an aggregate benchmark.', ...
    'instances', {{}}, 'error', '');
diary(fullfile(outdir, 'matlab.log'));
cleanup = onCleanup(@() diary('off'));
try
    assert(~isempty(which('yalmip')) && ~isempty(which('mosekopt')), 'YALMIP/MOSEK unavailable');
    addpath(fileparts(which('yalmip')), '-begin');
    frdir = fullfile(here, '..', 'third_party', 'frlib');
    addpath(frdir, fullfile(frdir, 'util'));
    assert(~isempty(which('frlibPrg')), 'frlib unavailable');
    ops0 = sdpsettings('solver','mosek','verbose',0, ...
        'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',120,'mosek.MSK_IPAR_NUM_THREADS',4);
    opsfr = sdpsettings('solver','frlib','frlib.solver','mosek', ...
        'frlib.approximation','dd','frlib.reduce','auto','verbose',0, ...
        'dualize',0, ...
        'savesolveroutput',1,'savesolverinput',1, ...
        'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',120,'mosek.MSK_IPAR_NUM_THREADS',4);
    names = {'N1_m5_s01', 'N1_m7_s01'};
    for ni = 1:numel(names)
        yalmip('clear');
        name = names{ni};
        fprintf('Instance %s\n', name);
        FD = jsondecode(fileread(fullfile(here, '..', 'data', 'benchmark', [name '.json'])));
        P = from_json(FD); r = str2double(FD.r); D = 4;
        entry = struct('name',name,'r',r,'D',D,'m',P.m,'q',floor(P.m/2), ...
            'basis_extraction_problem',NaN,'basis_extraction_info','', ...
            'gram_sizes',[],'candidate_monomials',[],'runs',{{}},'error','');
        try
            [~,~,info] = gen_sos_cert(P,D,r,ops0);
            entry.basis_extraction_problem = info.sol.problem;
            entry.basis_extraction_info = info.sol.info;
            assert(numel(info.v) == 3, 'Expected W/decrease/localizer Gram bases');
            vars = [info.x; info.y];
            E = info.E;
            mons = exponent_monomials(vars,E);
            bases = info.v;
            entry.gram_sizes = cellfun(@numel,bases);
            entry.candidate_monomials = E;
            entry.gram_block_order = {'W','decrease','localizer'};
            entry.gram_variable_order = [{'x'}, arrayfun(@(j) sprintf('y%d',j),1:P.N,'uni',0)];
            % Retain basis order and applied diagonals so the scaling can be traced.
            entry.gram_basis_exponents = cell(1,3);
            entry.restricted_jet_coefficients_exact = cell(1,3);
            Cs = cell(1,3); scales = cell(1,3); Escales = cell(1,3);
            for j = 1:3
                Ej = exponents_of(bases{j},vars);
                entry.gram_basis_exponents{j} = Ej;
                Cs{j} = restricted_jet(Ej,P.HS,entry.q);
                entry.restricted_jet_coefficients_exact{j} = symbolic_strings(Cs{j});
                [scales{j},Escales{j}] = jet_scales(Cs{j});
            end
            entry.jet_scale_targets_exact = Escales;
            entry.jet_scale_targets_binary64 = scales;
            % Rescale Gram bases while keeping the polynomial candidate spaces fixed.
            for scaled = [false true]
                c = sdpvar(numel(mons),1);
                W = c'*mons;
                dW = jacobian(W,vars)*info.f;
                Q = cell(1,3); sigma = cell(1,3); used = cell(1,3); constraints = [];
                for j = 1:3
                    if scaled, used{j}=scales{j}; else, used{j}=ones(numel(bases{j}),1); end
                    vj = used{j}.*bases{j};
                    Q{j} = sdpvar(numel(vj));
                    sigma{j} = vj'*Q{j}*vj;
                    constraints = [constraints, Q{j} >= 0];
                end
                k = min(2*P.m-2,D+P.m-1);
                cW = coefficients(W-sigma{1},vars);
                cD = coefficients(-dW-1e-3*info.x^k-sigma{2}-(r^2-info.x^2)*sigma{3},vars);
                constraints = [constraints, cW == 0, cD == 0];
                started = tic;
                sol = optimize(constraints,[],opsfr);
                elapsed = toc(started);
                if scaled, label='jet_row_scaled_dd'; else, label='unscaled_dd'; end
                row = struct('variant',label,'solver_problem',sol.problem,'solver_info',sol.info, ...
                    'wall_seconds',elapsed,'solver_seconds',sol.solvertime, ...
                    'reduced_gram_sizes',[],'coefficient_residual',NaN, ...
                    'classification','not classified','r_bad',NaN, ...
                    'witness',struct('x',NaN,'Y',[]), ...
                    'candidate_coefficients',[],'kernel_residual',[],'raw_file','');
                row.applied_scales_binary64 = used;
                row.prescribed_scales_exact = cell(1,3);
                for j = 1:3
                    if scaled, row.prescribed_scales_exact{j}=Escales{j};
                    else, row.prescribed_scales_exact{j}=repmat({'1'},size(used{j})); end
                end
                try
                    row.reduced_gram_sizes = sol.solveroutput.reducedmodel.K.s;
                catch
                end
                Qnum = cellfun(@value,Q,'uni',0);
                cval = value(c);
                if all(isfinite(cval))
                    row.candidate_coefficients = cval;
                    row.coefficient_residual = max(abs(value([cW;cD])));
                    kernels = zeros(1,2);
                    for j = 2:3
                        % Pull the Gram matrix back before measuring numerical kernel residuals.
                        Qoriginal = diag(used{j})*Qnum{j}*diag(used{j});
                        Cj = double(Cs{j});
                        kernels(j-1) = norm(Qoriginal*Cj,'fro')/max(1,norm(Qoriginal,'fro')*norm(Cj,'fro'));
                    end
                    row.kernel_residual = kernels;
                    Cinfo = struct('E',E,'cval',cval);
                    [st,rb,cand] = classify_candidate(sol.problem == 0,Cinfo,P,r,0);
                    row.classification=st; row.r_bad=rb;
                    row.witness=struct('x',cand.x,'Y',cand.Y);
                else
                    cand = struct;
                end
                raw = struct('name',name,'variant',label,'E',E,'c',cval, ...
                    'Q',{Qnum},'scales',{used},'r',r,'k',k,'classification',row.classification, ...
                    'witness',cand,'exact_scales',{Escales});
                raw.schema_version = 2;
                raw.gram_block_order = entry.gram_block_order;
                raw.gram_variable_order = entry.gram_variable_order;
                raw.gram_basis_exponents = entry.gram_basis_exponents;
                raw.restricted_jet_coefficients_exact = entry.restricted_jet_coefficients_exact;
                raw.jet_scale_targets_exact = Escales;
                raw.applied_scales_binary64 = used;
                raw.prescribed_scales_exact = row.prescribed_scales_exact;
                raw.exact_scales_role = ['Legacy field exact_scales contains the jet-derived symbolic targets ' ...
                    'for both variants; actual applied numerical diagonals are in scales/applied_scales_binary64.'];
                raw_path = fullfile(outdir,[name '_' label '.json']);
                write_json(raw_path,raw); row.raw_file=raw_path;
                entry.runs{end+1}=row;
                fprintf('%s: problem=%d, %s, %.2fs, reduced=%s\n',label,sol.problem, ...
                    row.classification,elapsed,mat2str(row.reduced_gram_sizes));
            end
        catch ex
            entry.error = getReport(ex,'extended','hyperlinks','off');
            fprintf('%s\n',entry.error);
        end
        report.instances{end+1}=entry;
        write_json(report_path,report);
    end
    if all(cellfun(@(z) isempty(z.error),report.instances))
        report.status='completed';
    else
        report.status='partial';
    end
catch ex
    report.status='blocked'; report.error=getReport(ex,'extended','hyperlinks','off');
    fprintf('%s\n',report.error);
end
write_json(report_path,report);
fprintf('Report: %s\n',report_path);
end

function P = from_json(F)
P=struct('N',F.N,'m',F.m);
for names={'Lam','B','a','v','w'}
    name=names{1}; val=read_exact(F.(name));
    P.([name 'S'])=val; P.(name)=double(val);
end
P.gS=read_exact(F.g); P.g=double(P.gS);
P.HS=read_exact(F.H_exact); P.H=double(read_exact(F.H_double));
end

function A = read_exact(F)
v=sym(zeros(numel(F.vals),1));
for i=1:numel(F.vals), v(i)=str2sym(F.vals{i}); end
A=reshape(v,double(F.dims(:)'));
end

function E = exponents_of(v,vars)
E=zeros(numel(v),numel(vars));
for i=1:numel(v)
    for j=1:numel(vars), E(i,j)=degree(v(i),vars(j)); end
end
end

function v = exponent_monomials(vars,E)
v=ones(size(E,1),1);
for j=1:numel(vars), v=v.*vars(j).^E(:,j); end
end

function C = restricted_jet(E,H,q)
xs=sym('xx'); C=sym(zeros(size(E,1),q+1));
h=H(:,1:q)*(xs.^((1:q)'));
for i=1:size(E,1)
    p=xs^E(i,1);
    for j=1:size(H,1), p=p*h(j)^E(i,j+1); end
    p=expand(p);
    for d=0:q, C(i,d+1)=subs(diff(p,xs,d),xs,0)/factorial(d); end
end
end

function [vals,exact] = jet_scales(C)
vals=ones(size(C,1),1); exact=cell(size(vals));
for i=1:size(C,1)
    [top,idx]=max(double(abs(C(i,:))));
    if top>0, scale=1/abs(C(i,idx)); else, scale=sym(1); end
    vals(i)=double(scale); exact{i}=char(scale);
end
end

function values = symbolic_strings(A)
values=cell(size(A));
for i=1:numel(A), values{i}=char(A(i)); end
end

function write_json(pathname,value)
fid=fopen(pathname,'w');
assert(fid>=0,'Could not open report output: %s',pathname);
closer=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value));
end
