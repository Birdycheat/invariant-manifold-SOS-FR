function report = run_reproduction(action, options)
%RUN_REPRODUCTION Reproduce the examples using MATLAB only.
arguments
    action (1,1) string = "quick"
    options.OutputDir (1,1) string = ""
    options.Benchmark (1,1) string = "all"
    options.InputDir (1,1) string = ""
    options.Instances = []
    options.Seeds (1,:) double = 1:20
    options.N (1,:) double = [1,3]
    options.M (1,:) double = [5,7]
    options.Enc (1,1) double = 1
    options.Round (1,1) logical = true
    options.KernelRatios (1,1) logical = false
    options.ScalarStrips (1,1) logical = false
end
root = fileparts(mfilename('fullpath'));
% Restore the caller's folder and path after every action.
oldPath = path; oldFolder = pwd;
cleanup = onCleanup(@() restoreSession(oldPath, oldFolder)); %#ok<NASGU>
addpath(genpath(fullfile(root, 'code')), genpath(fullfile(root, 'benchmarks')));
allowed = ["quick","symbolic","verify","anchors","witnesses","defects","tables","figures","solve","round","environment"];
assert(any(action == allowed), 'cm:Action', 'Action must be one of: %s.', strjoin(allowed, ', '));
environment = cm_environment(action);
outDir = cm_output_directory(root, options.OutputDir);
fprintf('Action: %s\nOutput: %s\n', action, outDir);
% Keep reference files unchanged while all outputs go to an external folder.
before = cm_manifest_check(root);
started = tic;
try
    switch action
        case "quick"
            result.symbolic = cm_symbolic_checks(root, outDir, 'smoke');
            result.certificates = cm_verify_certificates(root, outDir, 'smoke');
            result.tables = cm_tables(root, outDir);
        case "symbolic"
            result = cm_symbolic_checks(root, outDir, 'all');
        case "verify"
            result = cm_verify_certificates(root, outDir, 'all');
        case "anchors"
            result = cm_verify_anchors(root,outDir,'all');
        case "witnesses"
            result = cm_verify_witnesses(root,outDir,'all');
        case "defects"
            result = cm_verify_defects(root,outDir,'all');
        case "tables"
            result = cm_tables(root, outDir);
        case "figures"
            result = cm_readme_figures(root, outDir);
        case "solve"
            result = cm_solve_benchmark(root, outDir, char(options.Benchmark), options);
        case "round"
            assert(strlength(options.InputDir)>0, 'cm:InputDir', 'round requires InputDir.');
            result = cm_round_certificates(root, char(options.InputDir), outDir, char(options.Benchmark));
        case "environment"
            result = environment;
    end
    % Numerical completion and exact certification have distinct success criteria.
    if action=="quick"
        assert(strcmp(result.certificates.status,'passed'),'cm:Verification','A quick certificate check failed.');
        assert(result.symbolic.passed,'cm:Verification','A quick symbolic check failed.');
    elseif any(action==["verify","anchors","witnesses","defects","round"])
        assert(any(strcmpi(result.status,{'passed','PASS'})),'cm:Verification','The exact checks did not all pass. See the component report.');
    elseif action=="symbolic"
        assert(result.passed,'cm:Verification','An exact symbolic check failed. See the component report.');
    end
    after = cm_manifest_check(root);
    assert(isequal(before, after), 'cm:ReferenceChanged', 'The reference inventory changed.');
    status='passed';if action=="solve",status='completed';end
    report = struct('status',status,'action',char(action),'seconds',toc(started), ...
        'output_directory',outDir,'reference_unchanged',true,'environment',environment,'result',result);
catch err
    report = struct('status','failed','action',char(action),'seconds',toc(started), ...
        'output_directory',outDir,'identifier',err.identifier,'message',err.message);
    writeReport(outDir, report);
    rethrow(err);
end
writeReport(outDir, report);
fprintf('%s %s (%.2f s)\n', upper(report.status), action, report.seconds);
end

function writeReport(folder, report)
fid = fopen(fullfile(folder, 'run_report.json'), 'w', 'n', 'UTF-8');
assert(fid>=0, 'cm:Output', 'Cannot open the output report.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(report, PrettyPrint=true));
end

function restoreSession(oldPath, oldFolder)
path(oldPath);
if isfolder(oldFolder), cd(oldFolder); end
end
