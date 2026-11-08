function report = cm_verify_anchors(root,outDir,selection)
%CM_VERIFY_ANCHORS Prove the stored scalar anchored certificates on their strips.
if nargin<3,selection='all';end
outDir=cm_output_directory(root,outDir);
files=dir(fullfile(root,'data','benchmark','N*_m*_s*.json'));
if strcmp(selection,'smoke')
    wanted={'N1_m5_s01.json','N1_m7_s01.json','N3_m5_s01.json','N3_m7_s01.json'};
    files=files(ismember({files.name},wanted));
end
records=cell(numel(files),1);started=tic;
for k=1:numel(files)
    s=cm_read_json(fullfile(files(k).folder,files(k).name));clock=tic;
    rec=struct('instance',files(k).name,'proved',false,'failed_minor',0,'seconds',0);
    assert(strcmp(s.cm.status,'verifie'),'cm:Anchor','Stored anchored certificate is missing.');
    x=sym('x');ny=s.N;order=s.m;p=s.cm.p;
    Lam=cm_exact_flat(s.Lam);B=cm_exact_flat(s.B);a=cm_exact_flat(s.a);
    v=cm_exact_flat(s.v);w=cm_exact_flat(s.w);g=cm_exact_flat(s.g);g=g(:);
    H=cm_exact_flat(s.H_exact);mu=cm_exact(s.cm.mu);radius=cm_exact(s.cm.r);
    phi=H(:,1:p)*x.^((1:p).');dphi=diff(phi,x);
    F=expand(sum(g.*x.^((0:numel(g)-1).'))+x*a.'*phi);
    % Factor the prescribed central order and the transverse invariance defect.
    [qF,rF]=quorem(F,x^order,x);
    remainder=expand(dphi*F-Lam*phi-x*B*phi-v*x^2-w*x^3);
    qR=sym(zeros(ny,1));
    assert(isequal(rF,sym(0)),'cm:Anchor','The central factor is not divisible.');
    for j=1:ny
        [qR(j),rest]=quorem(remainder(j),x^(order-1),x);
        assert(isequal(rest,sym(0)),'cm:Anchor','The transverse defect is not divisible.');
    end
    M=Lam+x*B-x*dphi*a.';off=-(a-2*mu*qR)/2;
    G=[-qF,off.';off,-mu*(M+M.')];
    rec.proved=true;
    % A positive value at zero and no root prove each minor positive on the strip.
    for j=1:ny+1
        minor=expand(det(G(1:j,1:j)));
        if ~isAlways(subs(minor,x,0)>0) || cm_count_roots(minor,x,-radius,radius)>0
            rec.proved=false;rec.failed_minor=j;break;
        end
    end
    rec.seconds=toc(clock);records{k}=rec;
    fprintf('ANCHOR %d/%d %s: %d (%.2f s)\n',k,numel(files),files(k).name,rec.proved,rec.seconds);
end
count=sum(cellfun(@(r) r.proved,records));
report=struct('status','passed','selection',selection,'total',numel(files),'proved',count, ...
    'seconds',toc(started),'records',{records},'arithmetic','Exact rational Sylvester minors and closed-interval Sturm counts');
if count~=numel(files),report.status='failed';end
fid=fopen(fullfile(outDir,'anchors.json'),'w','n','UTF-8');guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
assert(count==numel(files),'cm:Anchor','An anchored certificate failed.');
end
