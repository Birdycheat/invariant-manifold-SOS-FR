function report = cm_verify_defects(root,outDir,selection)
%CM_VERIFY_DEFECTS Prove local invalidity of the scalar floating jet ablation.
% Reconstruct the invariant jet from the exact field. Candidate coefficients
% and its coordinate jet retain their exact IEEE binary64 interpretation.
if nargin<3,selection='all';end
selection=validatestring(selection,{'all','smoke'});
outDir=cm_output_directory(root,outDir);
files=dir(fullfile(root,'data','benchmark','N*_m*_s*.json'));
keep=false(size(files));
for k=1:numel(files),keep(k)=~isempty(regexp(files(k).name,'_m[57]_','once'));end
files=files(keep);
if strcmp(selection,'smoke')
    names={'N1_m5_s01.json','N1_m7_s01.json','N3_m5_s01.json','N3_m7_s01.json'};
    files=files(ismember({files.name},names));
end
methods={'jet sans RF','RF-jet'};records=cell(1,numel(files));sources=cell(size(records));
counts=zeros(2,4);started=tic;
for k=1:numel(files)
    source=fullfile(files(k).folder,files(k).name);S=cm_read_json(source);
    [H,premises]=reconstruct_jet(S);
    % Establish the finite hypotheses before applying the local invalidity criterion.
    assert(all(structfun(@(v) islogical(v)&&isscalar(v)&&v,premises)), ...
        'cm:DefectPremises','Finite invariant-manifold premises failed.');
    candidates=items(S.candidates);variants=cell(1,2);
    for method=1:2
        chosen=candidates(cellfun(@(c) strcmp(c.method,methods{method}),candidates));
        assert(numel(chosen)==1 && chosen{1}.solver_ok,'cm:DefectCandidate','Missing successful ablation output.');
        c=chosen{1};test=one_candidate(S,c,H);
        counts(method,1)=counts(method,1)+1;
        counts(method,2)=counts(method,2)+strcmp(c.status,'invalide');
        counts(method,3)=counts(method,3)+strcmp(c.status,'aucune violation');
        counts(method,4)=counts(method,4)+test.criterion_applies;
        variants{method}=struct('method',methods{method},'stored_status',c.status, ...
            'invalidity_established',test.criterion_applies,'transverse_defect',test);
    end
    records{k}=struct('instance',files(k).name,'m',S.m,'N',S.N, ...
        'jet_recomputed_from_field',true,'finite_premises',premises, ...
        'recomputed_h_through_m',{cellstr(string(H(:,2:end)))},'variants',{variants});
    sources{k}=struct('path',fullfile('data','benchmark',files(k).name),'sha256',cm_sha256(source));
    fprintf('DEFECTS %d/%d %s: %d/2 criteria established\n',k,numel(files),files(k).name, ...
        sum(cellfun(@(v) v.invalidity_established,variants)));
end
expected=80;if strcmp(selection,'smoke'),expected=4;end
assert(numel(files)==expected && all(counts(:,1)==expected),'cm:DefectCorpus','Unexpected scalar ablation corpus.');
summary=struct;
summary.without_FR=summary_row(counts(1,:));summary.with_FR=summary_row(counts(2,:));
if strcmp(selection,'all')
    assert(isequal(counts(:,2:3),[80 0;14 66]),'cm:DefectCorpus','Unexpected recorded witness counts.');
end
proved=sum(counts(:,4));
report=struct('status','passed','selection',selection,'instances',numel(files), ...
    'candidates',sum(counts(:,1)),'invalidity_established',proved,'summary',summary, ...
    'seconds',toc(started),'records',{records},'sources',{sources}, ...
    'source_sha256',cm_sha256([mfilename('fullpath') '.m']), ...
    'arithmetic','Exact rational field data; exact binary64 candidate and coordinate-jet values', ...
    'scope',['Applies the transverse-defect invalidity criterion after checking its finite premises. ' ...
        'The theorem is not re-proved here. No new witness search, PSD optimization or SDP is run. ' ...
        'rho values are numerical heuristic scales, not certified radii.']);
if proved~=sum(counts(:,1)),report.status='failed';end
fid=fopen(fullfile(outDir,'transversal_defects.json'),'w','n','UTF-8');assert(fid>=0);
closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
assert(strcmp(report.status,'passed'),'cm:DefectCriterion','A floating ablation output was not refuted by the criterion.');
end

function [H,checks]=reconstruct_jet(S)
% Coefficientwise invariance: Lambda h_d = -(known residual)_d.
n=double(S.N);m=double(S.m);A=cm_exact_flat(S.Lam);B=cm_exact_flat(S.B);
a=cm_exact_flat(S.a);v=cm_exact_flat(S.v);w=cm_exact_flat(S.w);
g=cm_exact_flat(S.g);g=g(:).';assert(numel(g)==m+1);
H=sym(zeros(n,m+1));
for degree=2:m
    [f,G,dH]=restricted_series(H,A,B,a,v,w,g,m);
    residual=G-dH*convolution_matrix(f,m);
    H(:,degree+1)=H(:,degree+1)-A\residual(:,degree+1);
end
[f,G,dH]=restricted_series(H,A,B,a,v,w,g,m);
residual=G-dH*convolution_matrix(f,m);
expected=sym(zeros(1,m+1));expected(end)=-1;
reference=cm_exact_flat(S.H_exact);
checks=struct('odd_m_at_least_three',mod(m,2)==1&&m>=3, ...
    'diagonal_Hurwitz_transverse_matrix',isequal(A,diag(diag(A)))&&all(isAlways(diag(A)<0)), ...
    'exact_reduced_leading_field_minus_x_power_m',isequal(f,expected), ...
    'finite_invariance_residual_through_m_zero',all(isAlways(residual(:)==0)), ...
    'h_has_no_constant_or_linear_term',all(isAlways(H(:,1:2)==0),'all'), ...
    'stored_jet_agrees_with_independent_reconstruction',isequal(H(:,2:end),reference(:,1:m)), ...
    'positive_prescribed_strip_radius',isAlways(cm_exact(S.r)>0));
end

function [f,G,dH]=restricted_series(H,A,B,a,v,w,g,m)
n=size(H,1);shifted=[sym(zeros(n,1)),H(:,1:m)];
f=g+a.'*shifted;G=A*H+B*shifted;G(:,3)=G(:,3)+v;G(:,4)=G(:,4)+w;
dH=[H(:,2:m+1)*diag(sym(1:m)),sym(zeros(n,1))];
end

function matrix=convolution_matrix(f,m)
matrix=sym(zeros(m+1));
for j=0:m,matrix(j+1,j+1:end)=f(1:m+1-j);end
end

function result=one_candidate(S,c,H)
n=double(S.N);m=double(S.m);p=double(c.p);q=(m-1)/2;limit=q+1;nv=n+1;
Hd=cm_exact(double(cm_exact_flat(S.H_double)),true);
assert(all(isAlways(Hd(:,1)==0)),'cm:CoordinateJet','The numerical coordinate jet must have no linear term.');
option=struct('H',Hd(:,1:p),'p',p);[fields,~,~]=cm_model('scalar',S,'xz',option);
field=cellfun(@(f) cm_poly('truncate',f,limit),fields{1},'UniformOutput',false);
E=double(cm_exact_flat(c.E));coefficients=cm_exact(double(cm_exact_flat(c.c)),true);coefficients=coefficients(:);
nonzero=~isAlways(coefficients==0);normalized=all(sum(E(nonzero,:),2)>=2);
W=cm_poly('truncate',cm_poly('from',E,coefficients),limit);D=cm_poly('zero',nv);
for j=1:nv
    product=cm_poly('mul',cm_poly('deriv',W,j),field{j});
    D=cm_poly('add',D,cm_poly('truncate',product,limit));
end
% Express the exact graph in the floating candidate's adapted coordinates.
replacement=cell(1,nv);replacement{1}=cm_poly('mono',1);
for j=1:n
    eta=H(j,2:q+1).';
    eta(1:min(p,q))=eta(1:min(p,q))-Hd(j,1:min(p,q)).';
    replacement{j+1}=cm_poly('from',(1:q).',eta);
end
% Extract the graph-restricted transverse gradient and transverse Hessian at zero.
delta=sym(zeros(n,q+1));N=sym(zeros(n));
for i=1:n
    restricted=cm_poly('compose',cm_poly('deriv',D,i+1),replacement);
    for degree=0:q,delta(i,degree+1)=coefficient_at(restricted,degree);end
    for j=1:n
        e=zeros(1,nv);e(i+1)=e(i+1)+1;e(j+1)=e(j+1)+1;
        N(i,j)=coefficient_at(D,e)*(1+(i==j));
    end
end
first=0;
for degree=1:q,if any(~isAlways(delta(:,degree+1)==0)),first=degree;break;end,end
checks=struct('normalized_candidate',normalized,'delta_0_zero',all(isAlways(delta(:,1)==0)), ...
    'N_negative_definite',cm_exact_pd(-N),'nonzero_low_order_delta',first>0);
applies=all(structfun(@(v) islogical(v)&&isscalar(v)&&v,checks));
infinity=zeros(1,q);euclidean=infinity;
for degree=1:q
    values=double(delta(:,degree+1));infinity(degree)=max(abs(values));euclidean(degree)=norm(values,2);
end
exponent=2./(m+1-2*(1:q));rhoInf=max(infinity.^exponent);rhoEuclid=max(euclidean.^exponent);
result=struct('checks',checks,'criterion_applies',applies,'first_nonzero_order',first, ...
    'N',{cellstr(string(N))},'delta_columns_0_through_q',{cellstr(string(delta))}, ...
    'heuristic_scale_infinity',rhoInf,'heuristic_scale_euclidean',rhoEuclid, ...
    'heuristic_scope','Numerical scales only; no certified witness radius is asserted.');
end

function value=coefficient_at(P,e)
[found,index]=ismember(e,P.E,'rows');value=sym(0);if found,value=P.c(index);end
end
function row=summary_row(counts)
row=struct('outputs',counts(1),'stored_explicit_witnesses',counts(2), ...
    'stored_no_witness_status',counts(3),'transverse_criterion_applies',counts(4), ...
    'invalidity_established',counts(4));
end
function rows=items(value)
if iscell(value),rows=value(:);else,rows=num2cell(value(:));end
end
