function report = cm_environment(action)
%CM_ENVIRONMENT Check the dependencies required by an action.
report = struct('MATLAB',version,'Symbolic',~isempty(which('sym')) && license('test','Symbolic_Toolbox')==1, ...
    'YALMIP',~isempty(which('sdpvar')),'MOSEK',~isempty(which('mosekopt')));
if any(string(action)==["quick","symbolic","verify","anchors","witnesses","defects","round"])
    assert(report.Symbolic, 'cm:SymbolicRequired', 'Install Symbolic Math Toolbox for this action.');
end
if string(action)=="solve"
    assert(report.Symbolic && report.YALMIP && report.MOSEK, 'cm:SolverRequired', ...
        'solve requires Symbolic Math Toolbox, YALMIP, MOSEK and a valid MOSEK license.');
end
end
