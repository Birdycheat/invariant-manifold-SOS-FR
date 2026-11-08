function report = cm_verify_system_classes(root, outDir, selection)
% Recompute the coupled/triple leading classes in exact rational arithmetic.
if nargin < 3, selection = 'all'; end
selection = validatestring(selection, {'smoke','all'});
records = {};
for bench = {'coupled','triple'}
    folder = fullfile(root, 'data', bench{1}, 'instances');
    files = dir(fullfile(folder, '*.json'));
    assert(~isempty(files), 'No %s instances found.', bench{1});
    [~, order] = sort({files.name}); files = files(order);
    if strcmp(selection, 'smoke'), files = files(1); end
    for k = 1:numel(files)
        source = jsondecode(fileread(fullfile(folder, files(k).name)));
        if strcmp(bench{1}, 'coupled'), row = planar(source);
        else, row = triple(source); end
        row.bench = bench{1}; row.instance = files(k).name;
        records{end+1} = row; %#ok<AGROW>
        fprintf('%s %s: %s.\n',bench{1},files(k).name,row.class);
    end
end
counts = struct('coupled_stable',0,'coupled_unstable',0, ...
    'triple_stable',0,'triple_unstable',0,'undecided',0);
for k = 1:numel(records)
    row = records{k}; field = [row.bench '_' row.class];
    if isfield(counts,field), counts.(field) = counts.(field) + 1;
    else, counts.undecided = counts.undecided + 1; end
end
report = struct('selection',selection, ...
    'scope',['Exact leading homogeneous classes, quartic manifold forcing, ' ...
             'and triangular Hurwitz transverse matrices; no SDP or enclosure test.'], ...
    'passed',~isempty(records) && all(cellfun(@(r) r.agrees,records)), ...
    'counts',counts,'records',{records});
if ~isempty(outDir)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    file = fullfile(outDir,'system_classes.json');
    fid = fopen(file,'w'); assert(fid >= 0, 'Cannot open report: %s',file);
    closer = onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(report));
end
fprintf('System classes: %d instances, passed=%d.\n',numel(records),report.passed);
end

function eta = source_eta(S, columns)
ny = double(S.ny);
lam = rational_matrix(S.Lam,ny,ny);
forcing = rational_matrix(S.Ehat,ny,columns);
% The quartic graph correction is recomputed before comparing stored values.
eta = -lam\forcing;
stored = rational_matrix(S.eta4,ny,columns);
assert(isequal(eta,stored), 'Stored eta4 differs from -Lam^{-1} Ehat.');
assert(all(all(tril(lam,-1) == 0)) && all(diag(lam) < 0), ...
    'Transverse matrix is not upper triangular with negative diagonal.');
end

function row = planar(S)
t = sym('t'); ny = double(S.ny);
eta = source_eta(S,5);
a = rational_matrix(S.a,2,ny);
g = rational_matrix(S.Fhat,2,6);
for j = 1:5
    for i = 1:2, g(i,j+i-1) = g(i,j+i-1) + a(i,:)*eta(:,j); end
end
assert(isequal(g,rational_matrix(S.lead,2,6)), ...
    'Stored leading field differs from the source equations.');
G = g*(t.^(0:5)).';
A = expand(G(1)+t*G(2)); B = expand(G(2)-t*G(1));
assert(~isequal(B,sym(0)), 'Every direction is invariant; this test does not decide it.');
signs = []; rays = {};
if coefficient(B,t,6) == 0
    radial = coefficient(A,t,6);
    assert(radial ~= 0, 'Degenerate vertical invariant ray.');
    sg = exact_sign(radial); signs(end+1) = sg;
    rays{end+1} = struct('vertical',true,'radial_sign',sg);
end
% Isolate angular roots exactly, then determine the radial sign on each ray.
[finite, iterations] = isolated_radial_signs(B,A,t);
for k = 1:numel(finite)
    signs(end+1) = finite{k}.radial_sign; %#ok<AGROW>
    rays{end+1} = finite{k}; %#ok<AGROW>
end
if isempty(signs)
    sequence = sturm_sequence(A,t); bound = root_bound(A,t);
    assert(variations(sequence,t,-bound)-variations(sequence,t,bound) == 0 ...
        && polynomialDegree(A,t) == 6 && coefficient(A,t,6) < 0 ...
        && subs(A,t,0) < 0, 'No-ray case not decided by strict radial negativity.');
    cls = 'stable'; method = 'Exact root exclusion proves globally negative radial polynomial';
elseif all(signs < 0)
    cls = 'stable'; method = 'Exact rational isolation of all invariant rays and strict radial signs';
else
    cls = 'unstable'; method = 'Exact rational isolation of an outward invariant ray';
end
row = struct('class',cls,'stored_class',S.leading.class,'method',method, ...
    'agrees',strcmp(cls,S.leading.class),'radial_polynomial',{coefficient_strings(A,t)}, ...
    'angular_polynomial',{coefficient_strings(B,t)},'rays',{rays}, ...
    'isolation_splits',iterations,'finite_distinct_rays',numel(finite));
end

function row = triple(S)
ny = double(S.ny); eta = source_eta(S,15);
a = rational_matrix(S.a,3,ny); g = rational_matrix(S.Fhat,3,21);
L4 = homogeneous_layer(4); L5 = homogeneous_layer(5); L8 = homogeneous_layer(8);
for i = 1:3
    for j = 1:size(L4,1)
        e = L4(j,:); e(i) = e(i)+1;
        [found,k] = ismember(e,L5,'rows'); assert(found);
        g(i,k) = g(i,k) + a(i,:)*eta(:,j);
    end
end
p = sym(zeros(size(L8,1),1));
for i = 1:3
    for j = 1:size(L5,1)
        e = L5(j,:); e(i) = e(i)+3;
        [found,k] = ismember(e,L8,'rows'); assert(found);
        p(k) = p(k)+g(i,j);
    end
end
pure = false(size(p)); diagonal = sym(zeros(1,3));
for i = 1:3
    e = zeros(1,3); e(i) = 8;
    [found,k] = ismember(e,L8,'rows'); assert(found);
    pure(k) = true; diagonal(i) = p(k);
end
% Weighted AM--GM bounds every mixed degree-eight term by pure powers.
weights = sym(zeros(1,3));
for i = 1:3
    for j = find(~pure).', weights(i) = weights(i)+abs(p(j))*sym(L8(j,i))/8; end
end
if all(diagonal < -weights)
    cls = 'stable';
elseif all(diagonal > weights)
    cls = 'unstable';
else
    cls = 'undecided';
end
row = struct('class',cls,'stored_class',S.leading.class, ...
    'method','Exact weighted AM-GM on the derivative of sum(x_i^4)/4', ...
    'diagonal',{symbolic_strings(diagonal)},'weights',{symbolic_strings(weights)}, ...
    'agrees',~strcmp(cls,'undecided') && strcmp(cls,S.leading.class));
end

function [rays, splits] = isolated_radial_signs(B,A,t)
% The Sturm sequence is computed once; every interval endpoint is a nonroot.
sequence = sturm_sequence(B,t); bound = root_bound(B,t);
left = variations(sequence,t,-bound); right = variations(sequence,t,bound);
queue = {struct('lo',-bound,'hi',bound,'vl',left,'vr',right)};
rays = {}; splits = 0;
while ~isempty(queue)
    box = queue{end}; queue(end) = [];
    count = box.vl-box.vr;
    if count == 0, continue; end
    assert(count > 0, 'Invalid Sturm variation count.');
    if count == 1
        [lo,hi] = polynomial_interval(A,t,box.lo,box.hi);
        if lo > 0 || hi < 0
            sg = 1; if hi < 0, sg = -1; end
            rays{end+1} = struct('vertical',false, ...
                'root_interval',{{char(box.lo),char(box.hi)}}, ...
                'radial_interval',{{char(lo),char(hi)}},'radial_sign',sg); %#ok<AGROW>
            continue
        end
    end
    divider = 2;
    mid = box.lo+(box.hi-box.lo)/divider;
    while subs(B,t,mid) == 0
        divider = divider+1; mid = box.lo+(box.hi-box.lo)/divider;
    end
    vm = variations(sequence,t,mid);
    queue{end+1} = struct('lo',box.lo,'hi',mid,'vl',box.vl,'vr',vm); %#ok<AGROW>
    queue{end+1} = struct('lo',mid,'hi',box.hi,'vl',vm,'vr',box.vr); %#ok<AGROW>
    splits = splits+1;
    assert(splits <= 10000, 'Radial sign could not be separated exactly.');
end
assert(numel(rays) == left-right, 'Not every distinct real ray was isolated.');
end

function sequence = sturm_sequence(p,x)
p = expand(p);
if polynomialDegree(p,x) == 0, sequence = {p}; return; end
[p,remainder] = quorem(p,gcd(p,diff(p,x)),x); assert(remainder == 0);
sequence = {expand(p),expand(diff(p,x))};
while true
    [~,remainder] = quorem(sequence{end-1},sequence{end},x);
    remainder = -expand(remainder);
    if remainder == 0, break; end
    sequence{end+1} = remainder; %#ok<AGROW>
end
end

function count = variations(sequence,x,at)
previous = 0; count = 0;
for k = 1:numel(sequence)
    sg = exact_sign(subs(sequence{k},x,at));
    if sg ~= 0
        if previous ~= 0 && sg ~= previous, count = count+1; end
        previous = sg;
    end
end
end

function bound = root_bound(p,x)
c = coeffs(p,x,'All'); bound = sym(0);
for k = 2:numel(c)
    candidate = abs(c(k)/c(1));
    if candidate > bound, bound = candidate; end
end
bound = 1+ceil(bound);
end

function [lo,hi] = polynomial_interval(p,x,left,right)
c = coeffs(p,x,'All'); lo = sym(0); hi = sym(0);
for k = 1:numel(c)
    products = [lo*left,lo*right,hi*left,hi*right];
    low = products(1); high = products(1);
    for j = 2:4
        if products(j) < low, low = products(j); end
        if products(j) > high, high = products(j); end
    end
    lo = low+c(k); hi = high+c(k);
end
end

function s = exact_sign(value)
if value > 0, s = 1; elseif value < 0, s = -1; else, s = 0; end
end

function c = coefficient(p,x,power)
c = subs(diff(p,x,power),x,0)/factorial(power);
end

function values = coefficient_strings(p,x)
values = symbolic_strings(coeffs(p,x,'All'));
end

function values = symbolic_strings(A)
values = arrayfun(@char,A,'UniformOutput',false);
end

function L = homogeneous_layer(degree)
L = zeros((degree+1)*(degree+2)/2,3); k = 0;
for a = degree:-1:0
    for b = degree-a:-1:0
        k = k+1; L(k,:) = [a,b,degree-a-b];
    end
end
end

function A = rational_matrix(value,rows,columns)
values = exact_leaves(value);
assert(numel(values) == rows*columns, 'Unexpected rational matrix dimensions.');
A = reshape(values,columns,rows).';
end

function values = exact_leaves(value)
if ischar(value) || (isstring(value) && isscalar(value))
    values = str2sym(char(value));
elseif iscell(value)
    values = sym(zeros(0,1));
    for row = 1:size(value,1)
        for column = 1:size(value,2)
            part = exact_leaves(value{row,column});
            values = [values;part(:)]; %#ok<AGROW>
        end
    end
elseif isnumeric(value)
    assert(all(isfinite(value(:))) && all(value(:) == fix(value(:))), ...
        'Expected exact integer or rational-string source coefficients.');
    values = sym(reshape(value.',[],1));
else
    error('Unsupported exact source encoding.');
end
end
