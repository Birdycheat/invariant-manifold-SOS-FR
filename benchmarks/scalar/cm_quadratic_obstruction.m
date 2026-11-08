function report = cm_quadratic_obstruction(root, outDir, selection)
% Necessary linear conditions for a quadratic candidate, not a feasibility test.
if nargin<3, selection='all'; end
files=dir(fullfile(root,'data','jet_inputs','*.json'));
if strcmp(selection,'smoke')
    names={'N1_m3_s01.json','N1_m5_s01.json','N3_m7_s01.json','ex1_D4.json'};
    files=files(ismember({files.name},names));
end
records=cell(1,numel(files)); counts=struct('infeasible_zero_only',0, ...
    'infeasible_indefinite_direction',0,'inconclusive_semidefinite_direction',0, ...
    'inconclusive_other',0);
for k=1:numel(files)
    S=cm_read_json(fullfile(files(k).folder,files(k).name));
    row=one_input(S); row.instance=files(k).name; records{k}=row;
    if isfield(counts,row.status), counts.(row.status)=counts.(row.status)+1;
    else, counts.inconclusive_other=counts.inconclusive_other+1; end
end
controls=psd_controls();
report=struct('count',numel(files),'counts',counts,'records',{records}, ...
    'psd_controls',{controls},'scope','Exact necessary conditions only; inconclusive does not imply feasible', ...
    'passed',~isempty(files));
if strcmp(selection,'all'), assert(report.count==124); else, assert(report.count==4); end
if ~isempty(outDir)
    fid=fopen(fullfile(outDir,'quadratic_obstruction.json'),'w'); assert(fid>=0);
    closer=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(report));
end
fprintf('Exact quadratic obstruction: %d inputs, %d zero-only and %d indefinite directions\n', ...
    report.count,counts.infeasible_zero_only,counts.infeasible_indefinite_direction);
end

function row = one_input(S)
ny=S.N; m=S.m; p=(m-1)/2; n=ny+1;
x=sym('x'); y=sym('y',[ny,1]); u=[x;y];
A=cm_exact_flat(S.Lam); B=cm_exact_flat(S.B); a=cm_exact_flat(S.a);
v=cm_exact_flat(S.v); w=cm_exact_flat(S.w); g=cm_exact_flat(S.g);
fx=g*(x.^((0:numel(g)-1).'))+x*(a.'*y);
fy=A*y+x*B*y+x^2*v+x^3*w;
h=sym(zeros(ny,1));
for degree=2:p
    residual=expand(subs(fy,y,h)-diff(h,x)*subs(fx,y,h));
    h=expand(h-(A\coefficient(residual,x,degree))*x^degree);
end
H=cm_exact_flat(S.H); expected=H(:,1:min(p,size(H,2)))*(x.^((1:min(p,size(H,2))).'));
zero(h-expected);
unknowns=sym('quadratic_coefficient',[n*(n+1)/2,1]); P=sym(zeros(n)); index=0;
for i=1:n
    for j=i:n, index=index+1; P(i,j)=unknowns(index); P(j,i)=unknowns(index); end
end
Wdot=2*u.'*P*[fx;fy]; restricted=expand(subs(jacobian(Wdot,y),y,h));
equations=sym(zeros(ny*(p+1),1)); index=0;
for i=1:ny
    for degree=0:p, index=index+1; equations(index)=coefficient(restricted(i),x,degree); end
end
% The nullspace gives necessary quadratic directions, not a feasibility proof.
constraints=jacobian(equations,unknowns); kernel=null(constraints);
direction={}; dimension=size(kernel,2);
if dimension==0
    status='infeasible_zero_only';
elseif dimension==1
    P0=subs(P,unknowns,kernel(:,1)); direction=arrayfun(@char,P0,'UniformOutput',false);
    if ~positive_semidefinite(P0) && ~positive_semidefinite(-P0)
        status='infeasible_indefinite_direction';
    else
        status='inconclusive_semidefinite_direction';
    end
else
    status=sprintf('inconclusive_dimension_%d',dimension);
end
row=struct('N',ny,'m',m,'jet_recomputed_from_field',true,'nullity',dimension, ...
    'status',status,'unique_direction',{direction});
end

function yes = positive_semidefinite(P)
assert(isequal(P,P.'));
yes=true;
% Semidefiniteness needs all principal minors, including nonleading ones.
for degree=1:size(P,1)
    ids=nchoosek(1:size(P,1),degree);
    for k=1:size(ids,1)
        if det(P(ids(k,:),ids(k,:)))<0, yes=false; return; end
    end
end
end

function records = psd_controls()
t=sym(1)/sym(10)^80;
matrices={sym([0 1;1 0]),[0 t;t 0],diag([t -t sym(0) sym(0)]), ...
    sym([1 1;1 1]),sym([-1 -1;-1 -1]),sym(zeros(4)),diag([t sym(1) sym(2) sym(3)]), ...
    [1 1;1 1-t],diag(sym([1 0 -1]))};
expected=[false false;false false;false false;true false;false true;true true;true false;false false;false false];
records=cell(1,numel(matrices));
for k=1:numel(matrices)
    actual=[positive_semidefinite(matrices{k}),positive_semidefinite(-matrices{k})];
    assert(isequal(actual,expected(k,:)));
    records{k}=struct('case',k,'psd',actual(1),'negative_psd',actual(2),'passed',true);
end
end

function c = coefficient(P,x,d)
c=sym(zeros(size(P)));
for k=1:numel(P)
    values=coeffs(expand(P(k)),x,'All'); if numel(values)>d, c(k)=values(end-d); end
end
end

function zero(P)
assert(all(isAlways(expand(P)==0),'all'),'Recomputed low jet differs from reference');
end
