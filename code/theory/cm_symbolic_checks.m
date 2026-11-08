function report = cm_symbolic_checks(root, outDir, selection)
% Exact symbolic reconstruction. Stored jets are comparison targets only.
if nargin < 3, selection = 'all'; end
selection = validatestring(selection, {'all','smoke'});
outDir = cm_output_directory(root, outDir);
report = struct('selection',selection,'matlab',version, ...
    'scope','Exact jets, exposing matrices, formal edge cases, BEG bounds and system classes; no SDP');
report.jets = cm_symbolic_jets(root, selection);
report.layers = cm_symbolic_layers(selection);
report.formal = cm_symbolic_formal();
report.quadratic = cm_quadratic_obstruction(root, outDir, selection);
report.beg = cm_symbolic_beg(root, selection);
report.classes = cm_verify_system_classes(root, outDir, selection);
report.passed = report.jets.passed && report.layers.passed && report.formal.passed && report.quadratic.passed && report.beg.passed && report.classes.passed;
fid = fopen(fullfile(outDir,'symbolic_checks.json'),'w');
assert(fid >= 0, 'Cannot open symbolic report');
closer = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
assert(report.passed, 'An exact symbolic check failed');
end
