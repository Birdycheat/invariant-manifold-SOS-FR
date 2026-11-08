function report = cm_solve_benchmark(root,outDir,bench,options)
%CM_SOLVE_BENCHMARK Execute numerical programs in a fresh external workspace.
% Options: Instances (file names), Seeds (default 1:20), N (default [1 3]),
% M (default [5 7]), Enc (default 1), Round (default true),
% KernelRatios/ScalarStrips (default false). Outputs never enter reference data.
arguments
    root (1,:) char
    outDir (1,:) char
    bench (1,:) char = 'all'
    options (1,1) struct = struct
end
allowed={'all','scalar','coupled','triple','beg','kernel-ratios','scalar-strips','dd-scaling'};
assert(ismember(bench,allowed),'Unknown benchmark.');
assert_external(outDir,root);
assert(license('test','Symbolic_Toolbox')==1,'Symbolic Math Toolbox is required.');
assert(exist('sdpvar','file')==2 && exist('mosekopt','file')~=0,'Configure YALMIP and MOSEK first.');
if ~isfolder(outDir),mkdir(outDir);end
% Copy producers and inputs into an external workspace before any numerical run.
work=fullfile(outDir,'workspace');assert(~isfolder(work),'Use a fresh output directory.');mkdir(work);
sourceCode=fullfile(root,'code');sourceBench=fullfile(root,'benchmarks');
assert(isfolder(sourceBench),'Missing benchmark MATLAB sources.');
codeFiles=dir(fullfile(sourceCode,'**','*.m'));benchFiles=dir(fullfile(sourceBench,'**','*.m'));
sourceNames=[{codeFiles.name},{benchFiles.name}];
assert(numel(unique(sourceNames))==numel(sourceNames),'cm:SourceCollision','MATLAB source names must be unique across benchmark and shared folders.');
% Flat temporary folders preserve the validated producer-relative data paths.
stage_sources(codeFiles,fullfile(work,'code'));
stage_sources(benchFiles,fullfile(work,'benchmarks'));
copyfile(fullfile(root,'third_party'),fullfile(work,'third_party'));
for family={'coupled','triple','beg'}
    target=fullfile(work,'data',family{1});mkdir(target);
    copyfile(fullfile(root,'data',family{1},'instances'),fullfile(target,'instances'));
end
for folder={'benchmark','jet_inputs','certificates','figdata'},mkdir(fullfile(work,'data',folder{1}));end
mkdir(fullfile(work,'tmp'));
if strcmp(bench,'dd-scaling')
    for name={'N1_m5_s01.json','N1_m7_s01.json'}
        copyfile(fullfile(root,'data','benchmark',name{1}),fullfile(work,'data','benchmark',name{1}));
    end
end
if strcmp(bench,'scalar-strips')
    inputs=dir(fullfile(root,'data','jet_inputs','*.json'));
    for j=1:numel(inputs),copyfile(fullfile(inputs(j).folder,inputs(j).name),fullfile(work,'data','jet_inputs',inputs(j).name));end
end
% Restore the MATLAB session and temporary-file locations even if a run fails.
oldPath=path;oldDir=pwd;oldTemp=getenv('TEMP');oldTmp=getenv('TMP');
cleanup=onCleanup(@() restore_session(oldPath,oldDir,oldTemp,oldTmp));
addpath(fullfile(work,'code'),fullfile(work,'benchmarks'),'-begin');
setenv('TEMP',fullfile(work,'tmp'));setenv('TMP',fullfile(work,'tmp'));
selected={bench};if strcmp(bench,'all'),selected={'scalar','coupled','triple','beg'};end
if option(options,'KernelRatios',false) && ~ismember('kernel-ratios',selected),selected{end+1}='kernel-ratios';end
if option(options,'ScalarStrips',false) && ~ismember('scalar-strips',selected),selected{end+1}='scalar-strips';end
records={};started=tic;
for j=1:numel(selected)
    family=selected{j};fprintf('SOLVE %s\n',family);clock=tic;
    folder=fullfile(outDir,family);if ~isfolder(folder),mkdir(folder);end
    execute_family(work,folder,family,options);
    rec=struct('benchmark',family,'status','completed','seconds',toc(clock));
    % Numerical completion and exact reconstruction have separate outcomes.
    if option(options,'Round',true) && ismember(family,{'scalar','coupled','triple','beg'})
        switch family
            case 'scalar',inputs=fullfile(work,'data','jet_inputs');
            case 'coupled',inputs=fullfile(folder,'center');
            case 'triple',inputs=fullfile(folder,'candidates');
            otherwise,inputs=folder;
        end
        rec.rounding=cm_round_certificates(work,inputs,fullfile(folder,'certificates'),family);
    end
    records{end+1}=rec; %#ok<AGROW>
end
report=struct('status','completed','benchmark',bench,'options',options,'workspace',work, ...
    'seconds',toc(started),'records',{records},'runtime','MATLAB only', ...
    'scope','New numerical points and timings can differ from the historical results.');
fid=fopen(fullfile(outDir,'solve_report.json'),'w','n','UTF-8');assert(fid~=-1);
done=onCleanup(@() fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end

function stage_sources(files,target)
mkdir(target);
for k=1:numel(files)
    copyfile(fullfile(files(k).folder,files(k).name),fullfile(target,files(k).name));
end
end

function execute_family(work,outputFolder,family,options)
cd(fullfile(work,'benchmarks'));
outdir=outputFolder;seeds=option(options,'Seeds',1:20); %#ok<NASGU>
seeds=seeds(:).';
assert(isnumeric(seeds)&&all(ismember(seeds,1:20)),'Seeds must be integers from 1 through 20.');
switch family
    case {'scalar','kernel-ratios'}
        Ns=option(options,'N',[1 3]);Ms=option(options,'M',[5 7]);
        Ns=Ns(:).';Ms=Ms(:).';
        assert(all(ismember(Ns,[1 3]))&&all(ismember(Ms,[3 5 7])),'Unsupported scalar N/M.');
        if strcmp(family,'kernel-ratios')
            assert(isequal(sort(Ns),[1 3])&&all(ismember([5 7],Ms))&&isequal(seeds,1:20), ...
                'Kernel-ratio regeneration requires all N=1,3, m=5,7 and seeds 1:20.');
        end
        cfgs=zeros(0,2);
        for Nsel=Ns
            for msel=Ms
                cfgs(end+1,:)=[Nsel msel]; %#ok<AGROW>
                run(fullfile(work,'benchmarks','run_batch_exact.m'));
            end
        end
        logf='jet_center_log.txt'; %#ok<NASGU>
        run(fullfile(work,'benchmarks','run_jet_exact.m'));
        if strcmp(family,'kernel-ratios')
            run(fullfile(work,'benchmarks','run_example1.m'));
            run(fullfile(work,'benchmarks','make_fig_data.m'));
        end
    case {'coupled','triple','beg'}
        instdir=fullfile(work,'data',family,'instances');
        files=selected_files(instdir,family,options);
        assert(~isempty(files),'No instances selected.');
        enc=option(options,'Enc',1); %#ok<NASGU>
        switch family
            case 'coupled',run(fullfile(work,'benchmarks','run_coupled_batch.m'));
            case 'triple',run(fullfile(work,'benchmarks','run_triple_batch.m'));
            case 'beg',run(fullfile(work,'benchmarks','run_beg_batch.m'));
        end
    case 'scalar-strips'
        Ns=option(options,'N',[1 3]);Ms=option(options,'M',[5 7]);
        [a,b]=ndgrid(Ns,Ms);cfgs=[a(:),b(:)];tag='selected'; %#ok<NASGU>
        run(fullfile(work,'benchmarks','run_scalar_regions.m'));
    case 'dd-scaling'
        review_dd_scaling(outputFolder,fullfile(outputFolder,'summary.json'));
end
diary off;
end

function files=selected_files(instdir,family,options)
L=dir(fullfile(instdir,'*.json'));files={L.name};requested=option(options,'Instances',{});
if ~isempty(requested)
    requested=cellstr(string(requested));
    for j=1:numel(requested),if ~endsWith(requested{j},'.json'),requested{j}=[requested{j} '.json'];end,end
    assert(all(ismember(requested,files)),'An explicitly selected instance is missing.');files=requested;
elseif ~strcmp(family,'beg') && isfield_case(options,'Seeds')
    seeds=option(options,'Seeds',1:20);keep=false(size(files));
    for j=1:numel(files),token=regexp(files{j},'_s(\d+)\.json$','tokens','once');keep(j)=~isempty(token)&&ismember(str2double(token{1}),seeds);end
    files=files(keep);
end
end
function value=option(s,name,fallback)
names=fieldnames(s);idx=find(strcmpi(names,name),1);if isempty(idx),value=fallback;else,value=s.(names{idx});end
end
function yes=isfield_case(s,name)
yes=any(strcmpi(fieldnames(s),name));
end
function assert_external(folder,root)
f=java.io.File(folder);p=char(f.getCanonicalPath());r=char(java.io.File(root).getCanonicalPath());
assert(f.isAbsolute(),'Output folder must be absolute.');
assert(~any(startsWith(regexp(lower(p),'[\\/]+','split'),'onedrive')),'Outputs inside OneDrive are forbidden.');
assert(~strcmpi(p,r)&&~startsWith(lower(p),[lower(r) filesep]),'Outputs must be outside the reference root.');
end
function restore_session(oldPath,oldDir,oldTemp,oldTmp)
path(oldPath);cd(oldDir);setenv('TEMP',oldTemp);setenv('TMP',oldTmp);diary off;
end
