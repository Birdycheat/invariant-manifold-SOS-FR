function report = cm_round_certificates(root,inputDir,outDir,bench)
%CM_ROUND_CERTIFICATES Rational rounding and exact affine/Gram projection.
% InputDir can be a raw JSON file or a directory. No numerical solver is used.
arguments
    root (1,:) char
    inputDir (1,:) char
    outDir (1,:) char
    bench (1,:) char {mustBeMember(bench,{'scalar','coupled','triple','beg'})}
end
assert_external(outDir,root);
if ~isfolder(outDir), mkdir(outDir); end
if isfile(inputDir)
    paths = {inputDir};
else
    listing = [dir(fullfile(inputDir,'**','*.json'));dir(fullfile(inputDir,'**','*.json.gz'))];
    paths = arrayfun(@(p) fullfile(p.folder,p.name),listing,'UniformOutput',false);
end
records = {}; started = tic;
for k = 1:numel(paths)
    item = cm_read_json(paths{k});
    if isfield(item,'center'), raw = item.center; else, raw = item; end
    if ~isfield(raw,'Q1') || ~isfield(raw,'E1'), continue; end
    if strcmp(bench,'beg') && ~isfield(item,'program') && ...
            ~contains(paths{k},'center_q2'), continue; end
    clock = tic;
    rec = struct('input',paths{k},'proved',false,'output','','reason','');
    try
        [cert,details] = round_one(root,item,raw,paths{k},bench);
        rec.proved = details.proved; rec.details = details;
        if rec.proved
            [~,name] = fileparts(paths{k});
            name = regexprep(name,'\.json$|\.cert$','');
            if strcmp(bench,'beg')
                [~,instance] = fileparts(char(cert.instance_file));
                name = [instance '_' char(cert.program)];
            end
            rec.output = fullfile(outDir,[name '.cert.json']);
            write_json(rec.output,cert);
        end
    catch err
        rec.reason = err.message;
        rec.error_identifier = err.identifier;
        rec.error_stack = err.stack;
    end
    rec.seconds = toc(clock); records{end+1} = rec; %#ok<AGROW>
    fprintf('ROUND %s: %s (%.2f s)\n',paths{k},pass_word(rec.proved),rec.seconds);
end
assert(~isempty(records),'No raw Gram certificate inputs found.');
proved = sum(cellfun(@(r) r.proved,records));
report = struct('status',pass_word(proved==numel(records)),'benchmark',bench, ...
    'inputs',numel(records),'proved',proved,'records',{records}, ...
    'seconds',toc(started),'arithmetic','Exact symbolic rational arithmetic', ...
    'scope','A failed projection is not a proof of infeasibility.');
write_json(fullfile(outDir,'rounding_report.json'),report);
end

function [cert,details] = round_one(root,item,raw,path,bench)
E = cell(1,3);
for b = 1:3
    data = raw.(sprintf('E%d',b));
    if isstruct(data), E{b} = double(cm_exact_flat(data));
    else, E{b} = double(cm_exact(data)); end
end
nv = size(E{1},2);
if strcmp(bench,'scalar')
    system = item; nx = 1; radius = cm_exact(item.r);
    [fields,~,~] = cm_model('scalar',system,'xz',struct('p',double(item.p)));
elseif strcmp(bench,'beg')
    nx = 1;
    [system,program,instanceFile] = beg_instance(root,item,path);
    if startsWith(program,'A_')
        radius = sym(1)/8; [fields,~,~] = cm_model('beg',system,'xz','surrogate');
    else
        denominator = str2double(extractAfter(program,'_r'));
        radius = sym(1)/denominator;
        entry = find_radius(system.radii,radius);
        [fields,~,~] = cm_model('beg',system,'xz',entry);
    end
else
    system = cm_read_json(fullfile(root,'data',bench,'instances',char(item.instance)));
    nx = 2 + strcmp(bench,'triple');
    if isfield(item,'r'), radius=cm_exact(item.r); else, radius=cm_exact(system.r); end
    if strcmp(bench,'coupled')
        [fields,~,~] = cm_model(bench,system,char(item.coords),double(item.enc));
    else
        [fields,~,~] = cm_model(bench,system,'xz');
    end
end
nvert = numel(fields);
gamma = cm_poly('mono',zeros(1,nv),radius^2);
domainVariables=nx;if strcmp(bench,'beg'),domainVariables=nv;end
for j = 1:domainVariables, e=zeros(1,nv);e(j)=2;gamma=cm_poly('add',gamma,cm_poly('mono',e),-1);end
if strcmp(bench,'scalar')
    e=zeros(1,nv);e(1)=double(item.k);ell=cm_poly('mono',e);
elseif strcmp(bench,'beg')
    if isfield(raw,'k'), degree=double(cm_exact(raw.k)); else, degree=8; end
    e=zeros(1,nv);e(1)=degree;ell=cm_poly('mono',e);
else
    ell=cm_poly('zero',nv);
    for j=1:nx,e=zeros(1,nv);e(j)=8;ell=cm_poly('add',ell,cm_poly('mono',e));end
end
if isfield(raw,'e'), eta=raw.e; else, eta=raw.eta; end
if strcmp(bench,'beg'), digits=[8 10 12 14]; else, digits=12; end
% Round first, then restore exact coefficient identities by rational projection.
for digit=digits
    Q1=rounded_matrix(raw.Q1,size(E{1},1),digit,strcmp(bench,'beg'));
    Q2=cell(1,nvert);Q3=Q2;
    for v=1:nvert
        Q2{v}=rounded_matrix(vertex_matrix(raw.Q2,v,nvert),size(E{2},1),digit,strcmp(bench,'beg'));
        Q3{v}=rounded_matrix(vertex_matrix(raw.Q3,v,nvert),size(E{3},1),digit,strcmp(bench,'beg'));
    end
    eta0=rounded_value(eta,digit,strcmp(bench,'beg'));
    if strcmp(bench,'beg')
        [Q1new,etaNew,faceResidual]=project_beg(E,Q1,eta0,fields,ell);
        W=cm_poly('gram',E{1},Q1new); S=cellfun(@(q) cm_poly('gram',E{3},q),Q3,'UniformOutput',false);
        Q3new=Q3;
    else
        composite=isfield(item,'composite') && logical(item.composite);
        [W,S,etaNew,faceResidual]=project_coefficients(E,Q1,Q3,eta0,fields,gamma,ell,nx,composite);
        Q1new=project_gram(E{1},Q1,W); Q3new=cell(1,nvert);
        for v=1:nvert,Q3new{v}=project_gram(E{3},Q3{v},S{v});end
    end
    Q2new=cell(1,nvert); identities=true;
    for v=1:nvert
        T=cm_poly('scale',lie(W,fields{v}),-1);
        T=cm_poly('add',T,ell,-etaNew);T=cm_poly('add',T,cm_poly('mul',S{v},gamma),-1);
        Q2new{v}=project_gram(E{2},Q2{v},T);
        identities=identities && cm_poly('equal',cm_poly('gram',E{2},Q2new{v}),T) ...
            && cm_poly('equal',cm_poly('gram',E{3},Q3new{v}),S{v});
    end
    identities=identities && cm_poly('equal',cm_poly('gram',E{1},Q1new),W);
    % An exact identity is certified only with positive Gram blocks and margin.
    matrices=[{Q1new},Q2new,Q3new]; pd=true;
    for j=1:numel(matrices),pd=pd && cm_exact_pd(matrices{j});end
    definite=pure_powers(E{1},1:nv) && pure_powers(E{2},nx+1:nv);
    ok=identities && pd && isAlways(etaNew>0) && definite;
    details=struct('proved',ok,'digits',digit,'exact_identities',identities, ...
        'positive_gram_blocks',pd,'positive_margin',isAlways(etaNew>0), ...
        'definite_bases',definite,'face_residual_before_projection',faceResidual, ...
        'vertices',nvert,'sizes',cellfun(@(a) size(a,1),E));
    if ok,break;end
end
if strcmp(bench,'scalar')
    keys={'N','m','p','D','k','r','Lam','B','a','v','w','g','H'};cert=struct;
    for j=1:numel(keys),cert.(keys{j})=item.(keys{j});end
    cert.Q1=text_matrix(Q1new);cert.Q2=text_matrix(Q2new{1});cert.Q3=text_matrix(Q3new{1});
else
    cert=struct('instance',item_instance(item,system),'r',char(radius));
    cert.Q1=triangle(Q1new);
    if strcmp(bench,'triple')
        cert.coordinates='x1, x2, x3, z = y - phi(x)';
        cert.Q2=triangle(Q2new{1});cert.Q3=triangle(Q3new{1});
    else
        cert.Q2=cellfun(@triangle,Q2new,'UniformOutput',false);
        cert.Q3=cellfun(@triangle,Q3new,'UniformOutput',false);
        if strcmp(bench,'beg'),cert.program=program;cert.instance_file=instanceFile;cert.k=degree_if_raw(raw);
        else,cert.coords=item.coords;cert.enc=item.enc;cert.composite=isfield(item,'composite')&&logical(item.composite);end
    end
end
for b=1:3,cert.(sprintf('E%d',b))=E{b};end
cert.eta=char(etaNew);cert.W=polynomial_rows(W);
end

function [W,S,eta,residual] = project_coefficients(E,Q1,Q3,eta,fields,gamma,ell,nx,composite)
WS=gram_support(E{1});SS=gram_support(E{3});MS=gram_support(E{2});
if composite,WS=WS(~(any(WS(:,1:nx),2)&any(WS(:,nx+1:end),2)),:);end
nW=size(WS,1);nS=size(SS,1);nv=size(E{1},2);nvert=numel(fields);n=nW+nvert*nS+1;
u=sym(zeros(n,1));u(1:nW)=coefficients_at(cm_poly('gram',E{1},Q1),WS);
for v=1:nvert,u(nW+(v-1)*nS+(1:nS))=coefficients_at(cm_poly('gram',E{3},Q3{v}),SS);end
u(end)=eta;rows=sym(zeros(0,n));
for v=1:nvert
    cols=cell(1,n);
    for j=1:nW,cols{j}=cm_poly('scale',lie(cm_poly('mono',WS(j,:)),fields{v}),-1);end
    for j=nW+1:n-1,cols{j}=cm_poly('zero',nv);end
    for j=1:nS,cols{nW+(v-1)*nS+j}=cm_poly('scale',cm_poly('mul',cm_poly('mono',SS(j,:)),gamma),-1);end
    cols{end}=cm_poly('scale',ell,-1);
    rows=[rows;outside_rows(cols,MS)]; %#ok<AGROW>
end
% Remove coefficients outside the permitted Gram support on every vertex.
[u,residual]=orthogonal_project(rows,u);
W=cm_poly('from',WS,u(1:nW));S=cell(1,nvert);
for v=1:nvert,S{v}=cm_poly('from',SS,u(nW+(v-1)*nS+(1:nS)));end
eta=u(end);
end

function [Q,eta,residual] = project_beg(E,Q,eta,fields,ell)
n=size(Q,1);[ii,jj]=find(triu(true(n)));idx=sortrows([ii,jj],[1 2]);
u=sym(zeros(size(idx,1)+1,1));cols=cell(1,numel(u));A=sym(zeros(0,numel(u)));
for j=1:size(idx,1),u(j)=Q(idx(j,1),idx(j,2));end
u(end)=eta;MS=gram_support(E{2});
for v=1:numel(fields)
    for j=1:size(idx,1)
        a=idx(j,1);b=idx(j,2);p=cm_poly('mono',E{1}(a,:)+E{1}(b,:),1+(a~=b));
        cols{j}=cm_poly('scale',lie(p,fields{v}),-1);
    end
    cols{end}=cm_poly('scale',ell,-1);A=[A;outside_rows(cols,MS)]; %#ok<AGROW>
end
[u,residual]=orthogonal_project(A,u);
for j=1:size(idx,1),Q(idx(j,1),idx(j,2))=u(j);Q(idx(j,2),idx(j,1))=u(j);end
eta=u(end);
end

function [u,residual] = orthogonal_project(A,u)
if isempty(A),residual=0;return;end
residual=max(abs(double(A*u)));
% Row reduction removes redundant constraints before the exact solve.
R=rref(A);R=R(any(R~=0,2),:);
if isempty(R),return;end
u=u-R.'*((R*R.')\(R*u));
assert(all(isAlways(A*u==0)),'Affine projection failed exactly.');
end

function A=outside_rows(cols,support)
allE=vertcat_cell_exponents(cols);bad=setdiff(unique(allE,'rows'),support,'rows');
A=sym(zeros(size(bad,1),numel(cols)));
for j=1:numel(cols),A(:,j)=coefficients_at(cols{j},bad);end
A=A(any(A~=0,2),:);
end

function allE=vertcat_cell_exponents(cols)
allE=zeros(0,size(cols{1}.E,2));
for j=1:numel(cols),allE=[allE;cols{j}.E];end %#ok<AGROW>
end

function values=coefficients_at(p,E)
values=sym(zeros(size(E,1),1));[found,at]=ismember(E,p.E,'rows');
if any(found),values(found)=p.c(at(found));end
end

function support=gram_support(E)
n=size(E,1);support=zeros(n*n,size(E,2));k=0;
for i=1:n,for j=1:n,k=k+1;support(k,:)=E(i,:)+E(j,:);end,end
support=unique(support,'rows');
end

function Q=project_gram(E,Q,target)
n=size(E,1);sums=zeros(n*n,size(E,2));k=0;
for j=1:n,for i=1:n,k=k+1;sums(k,:)=E(i,:)+E(j,:);end,end
[support,~,groups]=unique(sums,'rows');
assert(all(ismember(target.E,support,'rows')),'Target monomial outside the Gram space.');
t=coefficients_at(target,support);
for j=1:size(support,1),idx=find(groups==j);Q(idx)=Q(idx)+(t(j)-sum(Q(idx)))/numel(idx);end
end

function p=lie(W,field)
p=cm_poly('zero',size(W.E,2));
for j=1:numel(field)
    if any(W.E(:,j)>0),p=cm_poly('add',p,cm_poly('mul',cm_poly('deriv',W,j),field{j}));end
end
end

function Q=rounded_matrix(raw,n,d,beg)
if isstruct(raw)
    dimensions=double(raw.dims(:).');
    if isscalar(dimensions)
        dimensions=[dimensions 1];
    end
    Q=reshape(numeric_value(raw.vals),dimensions);
elseif iscell(raw) && numel(raw)==n && any(cellfun(@(q) numel(q)~=n,raw))
    Q=zeros(n);for i=1:n,Q(i,i:end)=numeric_value(raw{i});Q(i:end,i)=Q(i,i:end).';end
else
    Q=numeric_value(raw);
end
assert(isequal(size(Q),[n n]),'Gram matrix size mismatch.');
if beg,Q=(Q+Q.')/2;end
Q=rounded_value(Q,d,beg);Q=(Q+Q.')/2;
end

function q=rounded_value(value,d,beg)
if beg
    x=cm_exact(numeric_value(value),true)*sym(10)^d;
    lo=floor(x);frac=x-lo;y=lo+sym(isAlways(frac>sym(1)/2)) ...
        +sym(isAlways(frac==sym(1)/2)&isAlways(mod(lo,2)==1));
else
    x=numeric_value(value)*10^d;y=round(x);
    tie=abs(x-floor(x))==0.5;y(tie)=2*round(x(tie)/2);
    y=cm_exact(y,true);
end
q=y/sym(10)^d;
end

function value=numeric_value(raw)
% Numerical Gram input follows IEEE binary64 parsing before rational rounding.
if isnumeric(raw),value=double(raw);
elseif isa(raw,'sym'),value=double(raw);
elseif iscell(raw)
    if isvector(raw) && ~isempty(raw) && iscell(raw{1})
        raw=cm_json_rows(raw);
    end
    value=zeros(size(raw));
    for j=1:numel(raw),value(j)=numeric_value(raw{j});end
elseif isstring(raw) && ~isscalar(raw)
    value=zeros(size(raw));for j=1:numel(raw),value(j)=numeric_value(raw(j));end
elseif contains(char(raw),'/'),value=double(cm_exact(raw));
else,value=str2double(raw);
end
assert(all(isfinite(value(:))),'Non-finite numerical Gram input.');
end

function q=vertex_matrix(raw,v,nvert)
if nvert==1 && isstruct(raw),q=raw;return;end
if nvert==1 && ~iscell(raw)
    if ~ismatrix(raw),q=squeeze(raw(1,:,:));else,q=raw;end
    return;
end
if nvert==1 && iscell(raw) && numel(raw)>1,q=raw;return;end
if iscell(raw),q=raw{v};elseif isstruct(raw),q=raw(v);else,q=squeeze(raw(v,:,:));end
end

function ok=pure_powers(E,variables)
ok=true;
for j=variables,ok=ok&&any(E(:,j)>0 & sum(E,2)==E(:,j));end
end

function rows=text_matrix(Q)
rows=string(Q);
end
function rows=triangle(Q)
rows=cell(size(Q,1),1);for j=1:size(Q,1),rows{j}=string(Q(j,j:end));end
end
function rows=polynomial_rows(W)
rows=cell(size(W.E,1),1);for j=1:numel(rows),rows{j}={W.E(j,:),char(W.c(j))};end
end
function k=degree_if_raw(raw)
if isfield(raw,'k'),k=double(cm_exact(raw.k));else,k=8;end
end
function name=item_instance(item,system)
if isfield(item,'instance'),name=item.instance;else,name=system.a;end
end
function [I,program,instanceFile]=beg_instance(root,item,path)
if isfield(item,'program'),program=char(item.program);else,[~,program]=fileparts(path);end
L=dir(fullfile(root,'data','beg','instances','a_*.json'));I=[];
for j=1:numel(L)
    S=cm_read_json(fullfile(L(j).folder,L(j).name));
    if (isfield(item,'instance') && isAlways(cm_exact(item.instance)==cm_exact(S.a))) || contains(path,erase(L(j).name,'.json'))
        I=S;instanceFile=L(j).name;break;
    end
end
assert(~isempty(I),'Cannot identify the BEG parameter.');
end
function entry=find_radius(entries,r)
for j=1:numel(entries)
    if iscell(entries),e=entries{j};else,e=entries(j);end
    if isAlways(cm_exact(e.r)==r),entry=e;return;end
end
error('BEG radius not covered by an exact enclosure.');
end
function assert_external(path,root)
p=char(java.io.File(path).getCanonicalPath());r=char(java.io.File(root).getCanonicalPath());
assert(java.io.File(path).isAbsolute(),'Output folder must be absolute.');
assert(~any(startsWith(regexp(lower(p),'[\\/]+','split'),'onedrive')),'Outputs inside OneDrive are forbidden.');
assert(~strcmpi(p,r)&&~startsWith(lower(p),[lower(r) filesep]),'Outputs must be outside reference inputs.');
end
function word=pass_word(ok)
if ok,word='passed';else,word='failed';end
end
function write_json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid~=-1,'Cannot write report.');cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value));
end
