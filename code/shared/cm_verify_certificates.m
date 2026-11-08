function report = cm_verify_certificates(root, outDir, selection)
%CM_VERIFY_CERTIFICATES Verify retained certificates using exact MATLAB algebra.
if nargin < 3, selection = 'all'; end
root = char(java.io.File(char(root)).getCanonicalPath());
outDir = cm_output_directory(root,outDir);
allFiles = certificateFiles(root);
allRequested = (ischar(selection) || isstring(selection)) && isscalar(string(selection)) && strcmp(selection,'all');
if allRequested
    files = allFiles;
elseif (ischar(selection) || isstring(selection)) && isscalar(string(selection)) && strcmp(selection,'smoke')
    groups = {fullfile('data','certificates'),fullfile('data','certificates_regions'), ...
        fullfile('data','coupled','certificates'),fullfile('data','coupled','certificates_regions','jet4'), ...
        fullfile('data','triple','certificates')};
    files = cell(1,7);
    for k = 1:numel(groups)
        candidates = allFiles(startsWith(allFiles,[fullfile(root,groups{k}),filesep]));
        assert(~isempty(candidates)); files{k} = candidates{1};
    end
    beg = allFiles(contains(allFiles,[filesep,'beg',filesep,'certificates',filesep]));
    files{6} = beg{find(contains(beg,'_A_'),1)};
    files{7} = beg{find(contains(beg,'_B_'),1)};
else
    files = cellstr(string(selection));
    for k = 1:numel(files)
        if ~java.io.File(files{k}).isAbsolute(), files{k} = fullfile(root,files{k}); end
        assert(isfile(files{k}), 'cm:certificate:missing', 'Missing certificate: %s',files{k});
        files{k} = char(java.io.File(files{k}).getCanonicalPath());
    end
end
started = tic;
records = repmat(struct('file','','family','','valid',false,'checks',struct(), ...
    'seconds',0,'error',''),1,numel(files));
proofs = containers.Map('KeyType','char','ValueType','any');
% Recompute analytic bounds once per BEG parameter before checking its proofs.
if allRequested
    instances = dir(fullfile(root,'data','beg','instances','*.json'));
    for k = 1:numel(instances)
        S = cm_read_json(fullfile(instances(k).folder,instances(k).name));
        proof = cm_beg_bounds(S);
        assert(strcmp(proof.status,'passed') && proof.bounds_checked == 3);
        proofs(instances(k).name) = proof;
    end
end
for k = 1:numel(files)
    timer = tic;
    records(k).file = erase(files{k},[root,filesep]);
    try
        C = cm_read_json(files{k});
        if contains(files{k},[filesep,'coupled',filesep])
            family = 'coupled'; checks = coupledCertificate(root,C);
        elseif contains(files{k},[filesep,'triple',filesep])
            family = 'triple'; checks = tripleCertificate(root,C);
        elseif contains(files{k},[filesep,'beg',filesep])
            family = 'beg'; [checks,proofs] = begCertificate(root,files{k},C,proofs);
        else
            family = 'scalar'; checks = scalarCertificate(C);
        end
        records(k).family = family; records(k).checks = checks;
        records(k).valid = all(structfun(@(x) islogical(x) && isscalar(x) && x,checks));
    catch exception
        records(k).error = getReport(exception,'basic','hyperlinks','off');
    end
    records(k).seconds = toc(timer);
    if records(k).valid, status = 'PASS'; else, status = 'FAIL'; end
    fprintf('%s certificate %d/%d %s (%.2f s)\n',status,k,numel(files),records(k).file,records(k).seconds);
end
status = 'failed';
if all([records.valid]), status = 'passed'; end
bounds = 0; proofValues = values(proofs);
for k = 1:numel(proofValues), bounds = bounds+proofValues{k}.bounds_checked; end
report = struct('status',status,'total',numel(files),'passed',nnz([records.valid]), ...
    'analytic_beg_instances',proofs.Count,'analytic_beg_bounds',bounds, ...
    'arithmetic','exact symbolic rational coefficients and fraction-free principal minors', ...
    'sdp_solver_used',false,'seconds',toc(started),'certificates',records);
fid = fopen(fullfile(outDir,'certificates.json'),'w'); assert(fid >= 0);
guard = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); clear guard;
end

function files = certificateFiles(root)
groups = {fullfile('data','certificates','*.cert.json'),fullfile('data','certificates_regions','*.gz'), ...
    fullfile('data','coupled','certificates','*.gz'), ...
    fullfile('data','coupled','certificates_regions','jet4','*.gz'), ...
    fullfile('data','coupled','certificates_regions','composite','*.gz'), ...
    fullfile('data','triple','certificates','*.gz'),fullfile('data','beg','certificates','*.gz')};
files = {};
for k = 1:numel(groups)
    entries = dir(fullfile(root,groups{k}));
    files = [files,arrayfun(@(p) fullfile(p.folder,p.name),entries,'UniformOutput',false).']; %#ok<AGROW>
end
assert(numel(files)==376,'cm:certificate:inventory','Expected the 376 retained certificates.');
end

function checks = scalarCertificate(C)
nv = double(C.N)+1; E = {double(C.E1),double(C.E2),double(C.E3)};
Q = {cm_gram_matrix(C.Q1,'full'),cm_gram_matrix(C.Q2,'full'),cm_gram_matrix(C.Q3,'full')};
[fields,~,~] = cm_model('scalar',C,'xz');
r = cm_exact(C.r); eta = cm_exact(C.eta);
gamma = add(cm_poly('mono',zeros(1,nv),r^2),mono(nv,1,2),-1);
ell = mono(nv,1,double(C.k));
checks = baseChecks(E,Q,1,r,eta);
checks.margin_nonnegative = mod(double(C.k),2)==0 && C.k>0;
checks.identity = identities(E,Q(1),Q(2),Q(3),fields,gamma,ell,eta);
end

function checks = coupledCertificate(root,C)
S = cm_read_json(fullfile(root,'data','coupled','instances',C.instance));
nv = 2+double(S.ny); E = {double(C.E1),double(C.E2),double(C.E3)};
Q1 = cm_gram_matrix(C.Q1); Q2 = matrices(C.Q2); Q3 = matrices(C.Q3);
[fields,~,meta] = cm_model('coupled',S,C.coords,double(C.enc));
r = cm_exact(C.r); eta = cm_exact(C.eta);
gamma = cm_poly('mono',zeros(1,nv),r^2); ell = cm_poly('zero',nv);
for j = 1:2, gamma = add(gamma,mono(nv,j,2),-1); ell = add(ell,mono(nv,j,8)); end
checks = baseChecks(E,[{Q1},Q2,Q3],2,r,eta);
checks.four_vertices = numel(fields)==4 && numel(Q2)==4 && numel(Q3)==4;
checks.analytic_enclosures = any(meta.enclosure_order == [1,2]);
checks.identity = identities(E,{Q1},Q2,Q3,fields,gamma,ell,eta);
W = cm_poly('gram',E{1},Q1);
checks.recorded_polynomial = cm_poly('equal',W,storedPolynomial(C.W,nv));
if isfield(C,'composite') && C.composite
    checks.composite_structure = ~any(any(W.E(:,1:2)>0,2) & any(W.E(:,3:end)>0,2));
end
end

function checks = tripleCertificate(root,C)
S = cm_read_json(fullfile(root,'data','triple','instances',C.instance));
nv = 3+double(S.ny); E = {double(C.E1),double(C.E2),double(C.E3)};
Q = {cm_gram_matrix(C.Q1),cm_gram_matrix(C.Q2),cm_gram_matrix(C.Q3)};
[fields,~,~] = cm_model('triple',S,'xz'); r = cm_exact(S.r); eta = cm_exact(C.eta);
gamma = cm_poly('mono',zeros(1,nv),r^2); ell = cm_poly('zero',nv);
for j = 1:3, gamma = add(gamma,mono(nv,j,2),-1); ell = add(ell,mono(nv,j,8)); end
checks = baseChecks(E,Q,3,r,eta);
checks.identity = identities(E,Q(1),Q(2),Q(3),fields,gamma,ell,eta);
checks.recorded_polynomial = cm_poly('equal',cm_poly('gram',E{1},Q{1}),storedPolynomial(C.W,nv));
checks.stable_leading_class = strcmp(S.leading.class,'stable');
if isfield(C,'r'), checks.radius_matches = isAlways(cm_exact(C.r)==r); end
end

function [checks,proofs] = begCertificate(root,path,C,proofs)
[~,name,~] = fileparts(path); name = regexp(name,'^(.*?)_[AB]_','tokens','once');
instanceName = [name{1},'.json'];
S = cm_read_json(fullfile(root,'data','beg','instances',instanceName));
if ~isKey(proofs,instanceName)
    proof = cm_beg_bounds(S); assert(strcmp(proof.status,'passed') && proof.bounds_checked==3);
    proofs(instanceName) = proof;
end
proof = proofs(instanceName);
nv = 2; E = {double(C.E1),double(C.E2),double(C.E3)};
Q1 = cm_gram_matrix(C.Q1); Q2 = matrices(C.Q2); Q3 = matrices(C.Q3);
r = cm_exact(C.r); eta = cm_exact(C.eta); trueField = startsWith(C.program,'B_');
option = 'surrogate';
if trueField
    index = find(arrayfun(@(b) isAlways(cm_exact(b.r)==r),S.radii),1);
    assert(~isempty(index),'cm:certificate:radius','Unproved BEG radius.');
    option = S.radii(index);
end
[fields,~,~] = cm_model('beg',S,'xz',option);
gamma = add(add(cm_poly('mono',[0,0],r^2),mono(nv,1,2),-1),mono(nv,2,2),-1);
ell = mono(nv,1,double(C.k));
checks = baseChecks(E,[{Q1},Q2,Q3],1,r,eta);
checks.margin_nonnegative = mod(double(C.k),2)==0 && C.k>0;
checks.stable_class = strcmp(S.class,'stable');
checks.analytic_bounds = strcmp(proof.status,'passed') && proof.bounds_checked==3 ...
    && proof.T5_checked && proof.T11_checked;
checks.identity = identities(E,{Q1},Q2,Q3,fields,gamma,ell,eta);
if trueField
    checks.four_vertices = numel(fields)==4 && numel(Q2)==4 && numel(Q3)==4;
    % Corrupted coefficients, domains, margins and bounds must fail the identities.
    changed = fields;
    for k = 1:4, changed{k}{1} = add(changed{k}{1},cm_poly('mono',[1,1],sym(1)/1000)); end
    checks.control_coefficient = ~identities(E,{Q1},Q2,Q3,changed,gamma,ell,eta);
    larger = add(gamma,cm_poly('mono',[0,0],r^2*((sym(17)/16)^2-1)));
    checks.control_radius = ~identities(E,{Q1},Q2,Q3,fields,larger,ell,eta);
    checks.control_margin = ~identities(E,{Q1},Q2,Q3,fields,gamma,ell,2*eta);
    option.M = cm_exact(option.M)*sym(3)/4;
    changed = cm_model('beg',S,'xz',option);
    checks.control_bound = ~identities(E,{Q1},Q2,Q3,changed,gamma,ell,eta);
end
end

% Pure-power monomials control the SOS basis zero sets.
function checks = baseChecks(E,Q,nx,r,eta)
positive = true; symmetry = true;
for k = 1:numel(Q)
    symmetry = symmetry && isequal(Q{k},Q{k}.');
    positive = positive && cm_exact_pd(Q{k});
end
checks = struct('symmetry',symmetry,'positive_definite',positive, ...
    'epsilon_positive',isAlways(eta>0),'radius_positive',isAlways(r>0), ...
    'basis_zero_sets',purePowers(E{1},1:size(E{1},2)) && purePowers(E{2},nx+1:size(E{2},2)) ...
    && ~any(all(E{1}==0,2)));
end

function yes = purePowers(E,indices)
yes = true;
for j = indices, yes = yes && any(E(:,j)>0 & sum(E,2)==E(:,j)); end
end

function yes = identities(E,Q1,Q2,Q3,fields,gamma,ell,eta)
W = cm_poly('gram',E{1},Q1{1}); yes = numel(fields)==numel(Q2) && numel(Q2)==numel(Q3);
if ~yes, return; end
for k = 1:numel(fields)
    derivative = cm_poly('zero',size(E{1},2));
    for j = 1:numel(fields{k}), derivative = add(derivative,mul(cm_poly('deriv',W,j),fields{k}{j})); end
    % Compare every coefficient of -dot(W)-eta*ell-gamma*sigma with its SOS Gram.
    target = add(add(cm_poly('scale',derivative,-1),ell,-eta),mul(gamma,cm_poly('gram',E{3},Q3{k})),-1);
    yes = yes && cm_poly('equal',target,cm_poly('gram',E{2},Q2{k}));
end
end

function out = matrices(raw)
assert(iscell(raw)); out = cell(1,numel(raw));
for k = 1:numel(raw), out{k} = cm_gram_matrix(raw{k}); end
end

function p = storedPolynomial(raw,nv)
rows = cm_json_rows(raw); E = zeros(size(rows,1),nv); c = sym(zeros(size(rows,1),1));
for k = 1:size(rows,1), E(k,:) = double(rows{k,1}(:).'); c(k) = cm_exact(rows{k,2}); end
p = cm_poly('from',E,c);
end

function p = mono(nv,index,power)
E = zeros(1,nv); E(index) = power; p = cm_poly('mono',E);
end

function p = add(a,b,varargin)
p = cm_poly('add',a,b,varargin{:});
end

function p = mul(a,b)
p = cm_poly('mul',a,b);
end
