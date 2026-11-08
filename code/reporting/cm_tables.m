function result = cm_tables(root, outDir)
%CM_TABLES Reaggregate the retained numerical results; this does not prove PSD.
outDir=cm_output_directory(root,outDir);
methods = {'SOS D=2','SOS D=4','SOS D=6','Jones-Peet','frlib d','frlib dd','jet sans RF','RF-jet'};
labels = {'Free SOS D=2','Free SOS D=4','Free SOS D=6','Partially quadratic','Partial FR d','Partial FR dd','Jet without FR','Jet with FR d'};
data = dir(fullfile(root,'data','benchmark','N*_m*_s*.json'));
% The main scalar table excludes the additional order-three archive.
records = cell(1,numel(methods));
radii={};scalarMargins=[];scalarWitnesses=0;scalarExtraWitnesses=0;
for k=1:numel(data)
    s = cm_read_json(fullfile(data(k).folder,data(k).name));
    cands = list(s.candidates);
    found=sum(cellfun(@(c) strcmp(c.status,'invalide'),cands));scalarWitnesses=scalarWitnesses+found;
    if ~ismember(s.m,[5,7]),scalarExtraWitnesses=scalarExtraWitnesses+found;continue;end
    for j=1:numel(cands)
        c = cands{j}; idx = find(strcmp(methods,c.method));
        if isempty(idx), continue; end
        witness = strcmp(c.status,'invalide');
        ratio = NaN;
        if witness, ratio = 100*num(c.r_bad)/num(s.r); end
        records{idx}(end+1,:) = [s.N,s.m,logical(c.solver_ok),witness,strcmp(c.status,'aucune violation'),ratio];
    end
    [~,stem]=fileparts(data(k).name);
    certificate=cm_read_json(fullfile(root,'data','certificates',[stem '.cert.json']));
    [rounded,decreased]=radius_rounding(s.r,certificate.r);
    assert(rounded,'cm:Corpus','The certificate radius is not the downward 1e-4 rounding.');
    radii{end+1}=struct('instance',data(k).name,'input_radius',s.r, ...
        'certificate_radius',certificate.r,'exact_downward_grid_rounding',rounded, ...
        'decreased',decreased,'decrease',num(s.r)-num(certificate.r)); %#ok<AGROW>
    central=cm_read_json(fullfile(root,'data','jet_inputs',data(k).name));
    scalarMargins(end+1)=central.t; %#ok<AGROW>
end
scalarProtocol=struct('instances',numel(radii),'radii_decreased',sum(cellfun(@(r)r.decreased,radii)), ...
    'radii_unchanged',sum(cellfun(@(r)~r.decreased,radii)), ...
    'minimum_recorded_centering_margin',min(scalarMargins),'radius_records',{radii}, ...
    'radius_arithmetic','Exact unsigned-integer rational cross-products; no Symbolic toolbox', ...
    'margin_scope','Recorded floating optimization margin before rational reconstruction');
scalar = repmat(struct(),0,1);
for k=1:numel(methods)
    a = records{k}; group = [];
    for ny=[1,3]
        for order=[5,7]
            values = a(a(:,1)==ny & a(:,2)==order & a(:,4)==1,6);
            if ~isempty(values), group(end+1)=median(values); end %#ok<AGROW>
        end
    end
    row = struct('method',labels{k},'attempts',size(a,1),'solver_success',sum(a(:,3)), ...
        'witnesses',sum(a(:,4)),'none',sum(a(:,5)),'median_percent',med(a(a(:,4)==1,6)), ...
        'group_median_min_percent',min(group),'group_median_max_percent',max(group), ...
        'maximum_witness_percent',max(a(a(:,4)==1,6)));
    scalar=append(scalar,row);
end
assert(all([scalar.attempts]==80),'cm:Corpus','Expected 80 scalar inputs.');
multidimensional = struct(); centering = repmat(struct(),0,1);centeringProtocol=struct();costs=struct();
for bench={'coupled','triple'}
    name=bench{1}; s=cm_read_json(fullfile(root,'results',name,'classification.json')); rr=list(s.results);
    order={'SOS','frlib d','frlib dd','jet','Thm3'};
    if strcmp(name,'triple'), order={'SOS','frlib d','Thm2','jet','Thm3'}; end
    output=repmat(struct(),0,1);
    for k=1:numel(order)
        selected=rr(cellfun(@(r) strcmp(r.name,order{k}),rr));
        success=0; witnesses=0; none=0; ratios=[];
        for j=1:numel(selected)
            r=selected{j}; success=success+logical(r.ok);
            if strcmp(r.status,'invalid')
                witnesses=witnesses+1;
                input=cm_read_json(fullfile(root,'data',name,'candidates',r.instance));
                ratios(end+1)=100*num(r.r_bad)/num(input.r); %#ok<AGROW>
            elseif any(strcmp(r.status,{'no violation found','no_violation'})), none=none+1;
            end
        end
        row=struct('method',order{k},'attempts',numel(selected),'solver_success',success, ...
            'witnesses',witnesses,'none',none,'median_percent',med(ratios));
        output=append(output,row);
    end
    multidimensional.(name)=output;
    inputs=dir(fullfile(root,'data',name,'candidates','*.json')); before=[]; after=[]; xyBefore=[]; xyAfter=[];
    solverTimes=[];iterations=[];noDeletions=[];unstableMargins=[];
    for k=1:numel(inputs)
        s=cm_read_json(fullfile(inputs(k).folder,inputs(k).name));
        if ~strcmp(s.class,'stable'), unstableMargins(end+1)=s.t.xz_thm3;continue;end %#ok<AGROW>
        before(end+1)=s.t.xz_none; after(end+1)=s.t.xz_thm3; %#ok<AGROW>
        if isfield(s.t,'xy_thm2'), xyBefore(end+1)=s.t.xy_none;xyAfter(end+1)=s.t.xy_thm2; end %#ok<AGROW>
        if strcmp(name,'coupled'),meta=s.t;solverTimes(end+1)=s.t.center_time;
        else,meta=s.center;solverTimes(end+1)=s.center.time;end %#ok<AGROW>
        iterations(end+1)=meta.iter;noDeletions(end+1)=no_deletion(meta.deleted); %#ok<AGROW>
    end
    certificates=dir(fullfile(root,'data',name,'certificates','*.gz'));
    row=struct('benchmark',name,'coordinates','xz','stable_cases',numel(after), ...
        'maximum_abs_unreduced_margin',max(abs(before)),'minimum_reduced_margin',min(after), ...
        'exact_certificate_files',numel(certificates));
    centering=append(centering,row);
    centeringProtocol.(name)=struct('stable_cases',numel(iterations), ...
        'first_iteration_count',sum(iterations==1),'no_heuristic_deletion_count',sum(noDeletions), ...
        'all_first_iteration',all(iterations==1),'all_no_heuristic_deletion',all(noDeletions), ...
        'unstable_reduced_margins',statistics(unstableMargins), ...
        'scope','Retained integer metadata and floating optimization diagnostics; not a new solve');
    rounded=list(cm_read_json(fullfile(root,'results',name,'rounding.json')));
    roundingTimes=cellfun(@(r) r.seconds,rounded);
    certified=rounded(cellfun(@(r) r.proved,rounded));sizes=cellfun(@(r) r.sizes(:).',certified,'UniformOutput',false);
    sizes=vertcat(sizes{:});centeringProtocol.(name).maximum_recorded_rounded_gram_sizes=max(sizes,[],1);
    costs.(name)=struct('last_central_optimize_stable_seconds',statistics(solverTimes), ...
        'whole_rounding_all_instances_seconds',statistics(roundingTimes), ...
        'scope',['Central timer covers the final optimize call on stable inputs, excluding assembly and prior iterations. ' ...
        'Historical rounding timer covers the complete former per-instance routine on all retained inputs. ' ...
        'These are stored measurements, not timings of the native MATLAB implementation.']);
    if ~isempty(xyAfter)
        row=struct('benchmark',name,'coordinates','xy, ny=1','stable_cases',numel(xyAfter), ...
            'maximum_abs_unreduced_margin',max(abs(xyBefore)),'minimum_reduced_margin',min(xyAfter), ...
            'exact_certificate_files',NaN);
        centering=append(centering,row);
    end
end
s=cm_read_json(fullfile(root,'results','beg','classification.json')); rr=list(s.results);
programs={'A_SOS_D2','A_SOS_D4','A_SOS_D6','A_Jones_Peet','A_frlib_d','A_frlib_dd'};
selected=rr(cellfun(@(r) any(strcmp(r.program,programs)),rr));
beg=struct('attempts',numel(selected),'solver_success',sum(cellfun(@(r) logical(r.ok),selected)), ...
    'refuted',sum(cellfun(@(r) strcmp(r.status,'invalid'),selected)), ...
    'inconclusive',sum(cellfun(@(r) logical(r.ok) && ~strcmp(r.status,'invalid'),selected)));
assert(beg.refuted==78 && beg.inconclusive==1,'cm:Corpus','Unexpected BEG counts.');
begCentering=beg_centering(root);
% Count certificate problems separately from distinct system instances.
mainCohort=struct('scalar',numel(radii),'coupled',centering(1).exact_certificate_files, ...
    'triple',numel(dir(fullfile(root,'data','triple','certificates','*.gz'))), ...
    'beg_surrogate',begCentering.surrogate_certificate_files,'beg_true_parameters',begCentering.true_field_parameters);
mainCount=mainCohort.scalar+mainCohort.coupled+mainCohort.triple+mainCohort.beg_surrogate+mainCohort.beg_true_parameters;
standardScalar=sum([scalar(1:6).witnesses]);
standardCoupled=sum([multidimensional.coupled(1:3).witnesses]);
standardTriple=sum([multidimensional.triple(1:2).witnesses]);
aggregate=struct('standard_witness_outputs',standardScalar+standardCoupled+standardTriple+beg.refuted, ...
    'standard_breakdown',struct('scalar',standardScalar,'coupled',standardCoupled,'triple',standardTriple,'beg',beg.refuted), ...
    'main_exact_file_cohort',mainCohort,'main_exact_file_count',mainCount, ...
    'all_scalar_witness_outputs',scalarWitnesses,'additional_m3_scalar_witness_outputs',scalarExtraWitnesses, ...
    'scope','Output counts and file inventories; one largest retained BEG true-field radius per parameter, not distinct system count.');
result=struct('status','passed','scalar',scalar,'multidimensional',multidimensional,'centering',centering, ...
    'scalar_protocol',scalarProtocol,'centering_protocol',centeringProtocol,'recorded_costs',costs, ...
    'BEG_standard_search',beg,'BEG_centering',begCentering,'aggregate_counts',aggregate, ...
    'scalar_strip_comparison',scalar_regions(root),'gram_kernel_csv',gram_summary(root), ...
    'dd_scaling_evidence',cm_dd_evidence(root,outDir), ...
    'unit','100 * recorded r_bad / prescribed r, in percent', ...
    'scope','Stored-data aggregation. Feasibility, numerical centering and exact file inventories are separate.');
disp(struct2table(scalar));
if ~isfolder(outDir), mkdir(outDir); end
writetable(struct2table(scalar),fullfile(outDir,'table_I.csv'));
writetable(struct2table(multidimensional.coupled),fullfile(outDir,'table_II_coupled.csv'));
writetable(struct2table(multidimensional.triple),fullfile(outDir,'table_II_triple.csv'));
writetable(struct2table(centering),fullfile(outDir,'centering.csv'));
writetable(struct2table(vertcat(radii{:})),fullfile(outDir,'scalar_radius_rounding.csv'));
writetable(struct2table(result.gram_kernel_csv.groups),fullfile(outDir,'gram_kernel_groups.csv'));
writetable(struct2table(result.gram_kernel_csv.paper_pooled_groups),fullfile(outDir,'gram_kernel_pooled.csv'));
costRows=repmat(struct(),0,1);
for family={'coupled','triple'}
    c=costs.(family{1});s=c.last_central_optimize_stable_seconds;r=c.whole_rounding_all_instances_seconds;
    costRows=append(costRows,struct('benchmark',family{1},'stable_solver_cases',s.count, ...
        'median_last_optimize_seconds',s.median,'historical_rounding_cases',r.count, ...
        'median_historical_rounding_seconds',r.median));
end
writetable(struct2table(costRows),fullfile(outDir,'recorded_costs.csv'));
fid=fopen(fullfile(outDir,'tables.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function values=list(s)
if iscell(s),values=s(:);else,values=num2cell(s(:));end
end
function x=num(value)
if isnumeric(value),x=double(value);return;end
parts=split(string(value),'/');x=str2double(parts(1));
if numel(parts)==2,x=x/str2double(parts(2));end
end
function x=med(values)
if isempty(values),x=NaN;else,x=median(values);end
end
function rows=append(rows,row)
if isempty(rows),rows=row;else,rows(end+1)=row;end
end

function [rounded,decreased]=radius_rounding(before,after)
[bn,bd]=fraction(before);[an,ad]=fraction(after);
limit=intmax('uint64');grid=uint64(10000);
assert(bn<=idivide(limit,grid) && an<=idivide(limit,grid));
assert(an<=idivide(limit,bd) && bn<=idivide(limit,ad));
tick=idivide(bn*grid,bd,'floor');assert(tick<=idivide(limit,ad));
rounded=an*grid==tick*ad;decreased=an*bd<bn*ad;
end

function [n,d]=fraction(value)
value=char(string(value));
assert(~isempty(regexp(value,'^\d+(\.\d+)?(/\d+)?$','once')),'Unexpected nonnegative radius.');
parts=split(string(value),'/');numerator=char(parts(1));denominator=uint64(1);
dot=strfind(numerator,'.');
if ~isempty(dot),denominator=uint64(10^(numel(numerator)-dot));numerator=erase(numerator,'.');end
numeric=str2double(numerator);assert(numeric<=flintmax && numeric>=0 && numeric==fix(numeric));n=uint64(numeric);
if numel(parts)==2
    numeric=str2double(parts(2));assert(numeric>0 && numeric<=flintmax && numeric==fix(numeric));
    denominator=denominator*uint64(numeric);
end
d=denominator;
end

function yes=no_deletion(value)
yes=isempty(value) || (isnumeric(value) && all(value(:)==0));
end

function out=statistics(values)
values=values(:);values=values(isfinite(values));
if isempty(values),out=struct('count',0,'minimum',NaN,'median',NaN,'maximum',NaN);
else,out=struct('count',numel(values),'minimum',min(values),'median',median(values),'maximum',max(values));end
end

function out=scalar_regions(root)
lines=splitlines(string(fileread(fullfile(root,'results','coupled','regions_scalar.jsonl'))));
lines=lines(strlength(strtrim(lines))>0);ratios=[];floating=[];proved=0;
for k=1:numel(lines)
    row=jsondecode(lines(k));
    if row.proved,proved=proved+1;ratios(end+1)=row.r_proved/row.r_max_cm;end %#ok<AGROW>
    floating(end+1)=row.r_float/row.r_max_cm; %#ok<AGROW>
end
out=struct('attempts',numel(lines),'recorded_proved',proved, ...
    'proved_over_numerical_anchored_radius',statistics(ratios), ...
    'proved_near_ceiling_count',sum(ratios>=3.9),'floating_at_ceiling_count',sum(floating>=3.96), ...
    'scope',['Reaggregates retained bisection records; exact enlarged-strip certificates are checked separately. ' ...
    'The anchored reference is its numerical maximum, not an exactly certified maximum.'], ...
    'ceiling','Seven bisections on [r0,8*r0], anchored reference 2*r0; thresholds 3.9 for proved ratios and 3.96 for floating ratios');
end

function out=gram_summary(root)
table=readtable(fullfile(root,'data','figdata','gram_kernel.csv'));
keys=unique([table.N,table.m,table.i],'rows');groups=repmat(struct(),0,1);
for k=1:size(keys,1)
    keep=table.N==keys(k,1)&table.m==keys(k,2)&table.i==keys(k,3);
    groups=append(groups,struct('ny',keys(k,1),'m',keys(k,2),'jet_degree',keys(k,3), ...
        'count',nnz(keep),'median_ratio',median(table.rho(keep))));
end
keys=unique([table.m,table.i],'rows');pooled=repmat(struct(),0,1);forced=[];other=[];
for k=1:size(keys,1)
    keep=table.m==keys(k,1)&table.i==keys(k,2);value=median(table.rho(keep));
    pooled=append(pooled,struct('m',keys(k,1),'jet_degree',keys(k,2),'count',nnz(keep),'median_ratio',value));
    if keys(k,2)>=0 && keys(k,2)<=floor(keys(k,1)/2),forced(end+1)=value;
    elseif keys(k,2)==-1 || keys(k,2)==(keys(k,1)+1)/2,other(end+1)=value;end %#ok<AGROW>
end
out=struct('rows',height(table),'groups',groups,'paper_pooled_groups',pooled, ...
    'forced_group_medians',statistics(forced),'next_layer_and_random_group_medians',statistics(other), ...
    'scope','Stored CSV medians only. Unreduced historical Gram matrices are absent; this does not recompute Q*c norm ratios.');
end

function out=beg_centering(root)
points=cm_read_json(fullfile(root,'results','beg','central_points.json'));
files=dir(fullfile(root,'data','beg','certificates','*.gz'));iterations=[];noDeletions=[];margins=[];
best=containers.Map('KeyType','char','ValueType','double');surrogate=0;
for k=1:numel(files)
    parts=regexp(files(k).name,'^(.*?)_([AB]_.*?)\.cert\.json\.gz$','tokens','once');
    assert(~isempty(parts));meta=points.(parts{1}).(parts{2});
    iterations(end+1)=meta.iter;noDeletions(end+1)=no_deletion(meta.deleted);margins(end+1)=meta.t; %#ok<AGROW>
    if startsWith(parts{2},'A_'),surrogate=surrogate+1;
    else
        denominator=str2double(regexp(parts{2},'r(\d+)$','tokens','once'));
        radius=1/denominator;assert(ismember(denominator,[16 32 64]));
        if ~isKey(best,parts{1}) || radius>best(parts{1}),best(parts{1})=radius;end
    end
end
instances=dir(fullfile(root,'data','beg','instances','*.json'));stable=[];unstable=[];
for k=1:numel(instances)
    S=cm_read_json(fullfile(instances(k).folder,instances(k).name));[~,name]=fileparts(instances(k).name);
    value=points.(name).A_center_q2.t;
    if strcmp(S.class,'stable'),stable(end+1)=value;else,unstable(end+1)=value;end %#ok<AGROW>
end
names=sort(keys(best));radii=cellfun(@(name)best(name),names);chosen=cell(1,numel(names));
for k=1:numel(names),chosen{k}=struct('instance',names{k},'largest_radius',sprintf('1/%d',round(1/radii(k))));end
out=struct('retained_certificate_files',numel(files),'surrogate_certificate_files',surrogate, ...
    'true_field_certificate_files',numel(files)-surrogate,'true_field_parameters',numel(names), ...
    'largest_true_field_radius_counts',struct('r16',sum(radii==1/16),'r32',sum(radii==1/32),'r64',sum(radii==1/64)), ...
    'largest_true_field_radius_by_parameter',{chosen},'first_iteration_count',sum(iterations==1), ...
    'no_heuristic_deletion_count',sum(noDeletions),'all_first_iteration',all(iterations==1), ...
    'all_no_heuristic_deletion',all(noDeletions),'retained_centering_margins',statistics(margins), ...
    'surrogate_stable_centering_margins',statistics(stable),'surrogate_unstable_centering_margins',statistics(unstable), ...
    'scope','Recorded numerical center-point metadata before rational projection, not a new optimization.');
end
