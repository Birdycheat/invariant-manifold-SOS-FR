function report = cm_symbolic_beg(root, selection)
% Reconstruct BEG coefficients and analytic remainder bounds with rationals.
report.formulas = formulas();
files=dir(fullfile(root,'data','beg','instances','a_*.json'));
if strcmp(selection,'smoke')
    files=files(ismember({files.name},{'a_2_1.json','a_3_1.json','a_6_1.json'}));
end
records=cell(1,numel(files)); boundCount=0;
for k=1:numel(files)
    S=jsondecode(fileread(fullfile(files(k).folder,files(k).name)));
    record=cm_beg_bounds(S); record.instance=files(k).name;
    records{k}=record; boundCount=boundCount+record.bounds_checked;
    fprintf('Exact BEG reconstruction: %s, %d analytic bounds\n',files(k).name,record.bounds_checked);
end
report.instances=records; report.instance_count=numel(records); report.bound_count=boundCount; report.passed=true;
if strcmp(selection,'all'), assert(report.instance_count==15 && boundCount==45); end
if strcmp(selection,'smoke'), assert(report.instance_count==3 && boundCount==9); end
end

function out = formulas()
syms a positive
syms h2 h4
b=-a*(a-3)/(2*(a-1)); D0=1+1/(a-1);
u2=(b*h2+a^2/2)/a;
u4=(b*h4+b^2*h2^2/2+b*h2*a^2/2+a^4/24)/a;
% Solve the invariant-graph coefficients before extracting the reduced field.
Fy2=-h2+u2/D0;
h2v=simplify(solve(Fy2==0,h2)); assert(isscalar(h2v));
Fy4=subs(-h4-(u2^2-u4)/D0,h2,h2v);
h4v=simplify(solve(Fy4==0,h4)); assert(isscalar(h4v));
cubic=simplify(subs((a^3/6+a*b*h2-a*u2)/a,h2,h2v));
F5=(a^5/120+a^3*b*h2/6+a*(b*h4+b^2*h2^2/2)-u2*(a^3/6+a*b*h2)+a*(u2^2-u4))/a;
% The tricritical cubic term cancels; the quintic sign determines the listed classes.
quintic=simplify(subs(F5,[h2;h4],[h2v;h4v]));
rate=-1+b*(a-1)/a^2;
assert_zero(h2v-a/3); assert_zero(cubic);
assert_zero(quintic-a^4*(a^2-2*a-9)/(270*(a-1)^2));
assert_zero(rate+3*(a-1)/(2*a));
assert_zero(subs(Fy2,h2,h2v)); assert_zero(subs(Fy4,h4,h4v));
out=struct('domain','a > 1','b',char(b),'h2',char(h2v),'h4',char(h4v), ...
    'cubic',char(cubic),'c5',char(quintic),'transverse_eigenvalue',char(simplify(rate)), ...
    'threshold','1+sqrt(10); the order-five test is inconclusive at equality','passed',true);
end

function assert_zero(P)
assert(all(isAlways(simplify(expand(P))==0),'all'),'Exact polynomial comparison failed');
end

