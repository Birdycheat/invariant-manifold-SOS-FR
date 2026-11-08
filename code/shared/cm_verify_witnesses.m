function report = cm_verify_witnesses(root,outDir,selection)
%CM_VERIFY_WITNESSES Recheck stored invalidity witnesses in exact arithmetic.
if nargin<3,selection='all';end
outDir=cm_output_directory(root,outDir);records={};started=tic;
files=dir(fullfile(root,'data','benchmark','N*_m*_s*.json'));
if strcmp(selection,'smoke'),files=files(1);end
for k=1:numel(files)
    s=cm_read_json(fullfile(files(k).folder,files(k).name));
    [fields,~,~]=cm_model('scalar',s,'xy');candidates=items(s.candidates);
    % Reuse the candidate's binary64 coordinate jet, not the exact manifold jet.
    H=cm_exact(double(cm_exact_flat(s.H_double)),true);
    for j=1:numel(candidates)
        c=candidates{j};if ~strcmp(c.status,'invalide'),continue;end
        x=cm_exact(double(cm_exact(c.witness_x)),true);
        y=cm_exact(double(cm_exact_flat(c.witness_Y)),true);
        point=[x;y(:)];adapted=point;
        f=cellfun(@(p) cm_poly('eval',p,point),fields{1},'UniformOutput',false);f=[f{:}].';
        if c.p>0
            powers=sym(1:c.p);
            phi=H(:,1:c.p)*x.^powers.';
            dphi=H(:,1:c.p)*(powers.'.*x.^(powers.'-1));
            adapted(2:end)=y(:)-phi;f(2:end)=f(2:end)-dphi*f(1);
        end
        E=double(cm_exact_flat(c.E));coefficients=cm_exact(double(cm_exact_flat(c.c)),true);
        W=cm_poly('from',E,coefficients(:));gradient=gradientAt(W,adapted);
        lower=gradient.'*f;
        % Jones--Peet uses a full ball; the other scalar searches use a strip.
        if strcmp(c.method,'Jones-Peet')
            inside=isAlways(sum(point.^2)<=cm_exact(s.r)^2);
        else
            inside=isAlways(x^2<=cm_exact(s.r)^2);
        end
        records{end+1}=record('scalar',files(k).name,c.method,lower,inside); %#ok<AGROW>
    end
    fprintf('WITNESSES scalar %d/%d\n',k,numel(files));
end
for family={'coupled','triple'}
    bench=family{1};classification=cm_read_json(fullfile(root,'results',bench,'classification.json'));
    rows=items(classification.results);cache=containers.Map('KeyType','char','ValueType','any');
    models=containers.Map('KeyType','char','ValueType','any');
    if strcmp(selection,'smoke'),rows=rows(1);end
    for k=1:numel(rows)
        row=rows{k};if ~strcmp(row.status,'invalid'),continue;end
        key=char(row.instance);
        if ~isKey(cache,key)
            s=cm_read_json(fullfile(root,'data',bench,'instances',key));
            candidates=cm_read_json(fullfile(root,'data',bench,'candidates',key));
            cache(key)=struct('system',s,'candidates',{items(candidates.candidates)},'radius',candidates.r);
        end
        cached=cache(key);c=findCandidate(cached.candidates,row.name);
        E=double(c.E);coefficients=binaryCoefficients(c.c);W=cm_poly('from',E,coefficients);
        point=cm_exact(row.witness);point=point(:);
        modelKey=[key '#' c.coords];
        if ~isKey(models,modelKey)
            [fields,~,meta]=cm_model(bench,cached.system,c.coords,1);
            models(modelKey)=struct('fields',{fields},'meta',meta);
        end
        model=models(modelKey);meta=model.meta;
        if strcmp(bench,'coupled')
            gradient=gradientAt(W,point);parts=meta.parts;values=sym(zeros(1,3));
            for a=1:3
                f=cellfun(@(p) cm_poly('eval',p,point),parts{a},'UniformOutput',false);
                values(a)=gradient.'*[f{:}].';
            end
            lower=coupledLower(values,point(1)-point(2),sum(point(1:2).^2));nx=2;
        else
            lower=dotAt(W,model.fields{1},point);nx=meta.nx;
        end
        inside=isAlways(sum(point(1:nx).^2)<=cm_exact(cached.radius)^2);
        records{end+1}=record(bench,key,row.name,lower,inside); %#ok<AGROW>
    end
    fprintf('WITNESSES %s completed\n',bench);
end
classification=cm_read_json(fullfile(root,'results','beg','classification.json'));rows=items(classification.results);
if strcmp(selection,'smoke'),rows=rows(1);end
cache=containers.Map('KeyType','char','ValueType','any');
for k=1:numel(rows)
    row=rows{k};if ~strcmp(row.status,'invalid'),continue;end
    key=char(row.instance);
    if ~isKey(cache,key)
        s=cm_read_json(fullfile(root,'data','beg','instances',[key '.json']));
        candidates=cm_read_json(fullfile(root,'data','beg','candidates',[key '.json']));
        [proof,bounds]=cm_beg_bounds(s);assert(strcmp(proof.status,'passed'));
        cache(key)=struct('system',s,'candidates',{items(candidates.candidates)},'bounds',bounds);
    end
    cached=cache(key);s=cached.system;c=findCandidate(cached.candidates,row.program);
    W=cm_poly('from',double(c.E),binaryCoefficients(c.c));
    point=cm_exact(row.witness);point=point(:);
    if ~strcmp(c.coords,'xy'),point(2)=point(2)-cm_exact(s.h2)*point(1)^2;end
    if startsWith(row.program,'A_')
        [fields,~,~]=cm_model('beg',s,c.coords,'surrogate');
        lower=dotAt(W,fields{1},point);radius=sym(1)/8;
    else
        % Subtract the analytic remainder to obtain a rigorous derivative lower bound.
        gradient=gradientAt(W,point);fields=begRows(s.T11xz);
        lower=dotAt(W,fields,point)-sum(abs(gradient))*cached.bounds.bounds(1).B*sum(point.^2)^6;
        radius=cached.bounds.bounds(1).r;
    end
    records{end+1}=record('beg',key,row.program,lower,isAlways(sum(point.^2)<radius^2)); %#ok<AGROW>
end
dd=cm_read_json(fullfile(root,'results','generalisation','dd_scaling_review.json'));
ddInstances=items(dd.instances);
if strcmp(selection,'smoke'),ddInstances=ddInstances(1);end
for k=1:numel(ddInstances)
    instance=ddInstances{k};s=cm_read_json(fullfile(root,'data','benchmark',[instance.name '.json']));
    [fields,~,~]=cm_model('scalar',s,'xy');runs=items(instance.runs);
    for j=1:numel(runs)
        run=runs{j};W=cm_poly('from',double(instance.candidate_monomials),cm_exact(run.candidate_coefficients,true));
        point=cm_exact([run.witness.x;run.witness.Y(:)],true);
        lower=dotAt(W,fields{1},point);
        records{end+1}=record('dd-scaling',instance.name,run.variant,lower,isAlways(point(1)^2<=cm_exact(s.r)^2)); %#ok<AGROW>
    end
end
passed=sum(cellfun(@(r) r.proved,records));
report=struct('status','passed','selection',selection,'total',numel(records),'passed',passed, ...
    'seconds',toc(started),'records',{records}, ...
    'scope','Stored witnesses only. A missing witness is inconclusive; no new search is performed.');
if passed~=numel(records),report.status='failed';end
fid=fopen(fullfile(outDir,'witnesses.json'),'w','n','UTF-8');guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
assert(passed==numel(records),'cm:Witness','A stored witness did not pass its exact sign or domain check.');
end

function result=record(family,instance,method,lower,inside)
result=struct('family',family,'instance',instance,'method',method,'proved',inside && isAlways(lower>0), ...
    'lower_bound',char(lower),'inside_domain',inside);
end
function rows=items(value)
if iscell(value),rows=value(:);else,rows=num2cell(value(:));end
end
function c=findCandidate(candidates,name)
index=find(cellfun(@(c) strcmp(c.name,name),candidates));assert(isscalar(index));c=candidates{index};
end
function values=binaryCoefficients(c)
values=cm_exact(double(cm_exact(c)),true);values=values(:);
end
function gradient=gradientAt(W,point)
gradient=sym(zeros(numel(point),1));
for j=1:numel(point),gradient(j)=cm_poly('eval',cm_poly('deriv',W,j),point);end
end
function value=dotAt(W,field,point)
gradient=gradientAt(W,point);value=sym(0);
for j=1:numel(point),value=value+gradient(j)*cm_poly('eval',field{j},point);end
end
% Alternating Taylor intervals certify signs without floating transcendental calls.
function lower=coupledLower(values,u,rho)
lower=sym(0);
for order=[4,8,16,32,64]
    cosine=sym(1);term=sym(1);
    for k=1:order,term=-term*u^2/((2*k-1)*(2*k));cosine=cosine+term;end
    cosineError=abs(term*u^2)/((2*order+1)*(2*order+2));
    exponential=sym(1);term=sym(1);
    for k=1:order,term=-term*rho/k;exponential=exponential+term;end
    exponentialError=abs(term*rho)/(order+1);
    bounds=[2-cosine-cosineError,2-cosine+cosineError;exponential-exponentialError,exponential+exponentialError];
    lower=values(1);
    for j=1:2
        endpoint=1+isAlways(values(j+1)<0);lower=lower+values(j+1)*bounds(j,endpoint);
    end
    if isAlways(lower>0),return;end
end
end
function fields=begRows(value)
rows=cm_json_rows(value);fields=cell(1,2);
for j=1:2
    selected=cellfun(@(x) double(cm_exact(x))==j,rows(:,1));
    E=double(cm_exact(rows(selected,2:3)));c=cm_exact(rows(selected,4));
    fields{j}=cm_poly('from',E,c);
end
end
