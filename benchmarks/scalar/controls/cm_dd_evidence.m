function report = cm_dd_evidence(root,outDir)
% Check the retained DD control without Symbolic Math Toolbox or an SDP.
outDir = cm_output_directory(root,outDir);
arithmeticChecks = arithmetic_self_checks();
summaryPath = fullfile(root,'results','generalisation','dd_scaling_review.json');
indexPath = fullfile(root,'results','generalisation','dd_scaling_traceability_review2.json');
summary = cm_read_json(summaryPath); index = cm_read_json(indexPath);
assert(strcmp(summary.status,'completed'),'The retained DD experiment did not complete.');
check_hash(summaryPath,index.historical_summary);
rawFiles = items(index.raw_files); sidecars = items(index.sidecars);
raw = containers.Map('KeyType','char','ValueType','any');
% Authenticate the recovered raw records before comparing their contents.
for k = 1:numel(rawFiles)+numel(sidecars)
    if k <= numel(rawFiles), mapping = rawFiles{k};
    else, mapping = sidecars{k-numel(rawFiles)}; end
    file = contained_path(root,mapping.copy_path); check_hash(file,mapping);
    if k <= numel(rawFiles)
        key = [mapping.instance '#' mapping.variant];
        assert(~isKey(raw,key),'Duplicate DD raw-file mapping.');
        raw(key) = struct('data',cm_read_json(file),'mapping',mapping);
    end
end
assert(raw.Count == 4,'Expected four retained raw DD files.');
checks = {}; entries = items(summary.instances);
for k = 1:numel(entries)
    entry = entries{k}; assert(isempty(entry.error));
    source = cm_read_json(fullfile(root,'data','benchmark',[entry.name '.json']));
    runs = items(entry.runs); assert(numel(runs) == 2);
    for j = 1:numel(runs)
        run = runs{j};
        assert(run.solver_problem == 0 && strcmp(run.classification,'invalide'));
        % Recompute the witness derivative with exact rational arithmetic.
        [x,dW] = derivative(source,entry.candidate_monomials,run);
        assert(qcompare(qabs(x),qtext(source.r)) < 0,'DD witness lies outside the strip.');
        assert(dW.n.signum() > 0,'DD witness is not an exact positive derivative.');
        saved = raw([entry.name '#' run.variant]);
        coherence = raw_consistency(saved.data,saved.mapping,entry,run,source);
        checks{end+1} = struct('instance',entry.name,'variant',run.variant, ...
            'inside_strip',true,'dW_strictly_positive_exact',true, ...
            'dW_exact',qstring(dW),'raw_traceability',coherence); %#ok<AGROW>
    end
end
assert(numel(checks) == 4);
report = struct('status','passed', ...
    'arithmetic','Exact BigInteger rationals and direct binary64 decoding; no Symbolic toolbox', ...
    'checks_passed',numel(checks),'recovered_raw_files_checked',raw.Count, ...
    'exact_arithmetic_self_checks',arithmeticChecks, ...
    'historical_summary_unchanged',true,'traceability_index_sha256',cm_sha256(indexPath), ...
    'scope',['Exact invalidity witnesses and retained raw-file coherence. Historical Gram ' ...
             'basis order and reduced-model construction remain unverified because they were not retained.'], ...
    'checks',{checks});
file = fullfile(outDir,'dd_evidence.json'); fid = fopen(file,'w');
assert(fid >= 0,'Cannot open DD evidence report.'); closer = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
fprintf('DD evidence: four exact witnesses and four recovered raw files passed.\n');
end

function coherence = raw_consistency(raw,mapping,entry,run,source)
assert(strcmp(raw.name,entry.name) && strcmp(raw.variant,run.variant));
assert(isequal(raw.E,entry.candidate_monomials));
assert(isequal(vector(raw.c),vector(run.candidate_coefficients)));
assert(raw.r == entry.r && raw.k == min(2*entry.m-2,entry.D+entry.m-1));
assert(strcmp(raw.classification,run.classification));
assert(raw.witness.x == run.witness.x && isequal(vector(raw.witness.Y),vector(run.witness.Y)));
assert(strcmp(mapping.source_path,run.raw_file));
blocks = items(raw.Q); scales = items(raw.scales); targets = items(raw.exact_scales);
sizes = zeros(1,numel(blocks)); applied = cell(1,numel(blocks));
for k = 1:numel(blocks)
    Q = double(blocks{k}); sizes(k) = size(Q,1);
    assert(isequal(size(Q),[sizes(k),sizes(k)]) && all(isfinite(Q(:))));
    assert(isequal(Q,Q.'),'A retained DD Gram matrix is not exactly symmetric.');
    values = vector(scales{k}); exact = text_vector(targets{k});
    assert(numel(values) == sizes(k) && numel(exact) == sizes(k));
    assert(all(isfinite(values)) && all(values > 0),'Applied DD scales must be positive and finite.');
    for i = 1:sizes(k)
        if strcmp(run.variant,'unscaled_dd'), prescribed = qinteger(1);
        else, prescribed = qtext(exact{i}); end
        assert(rounds_to(prescribed,values(i)), ...
            'Applied DD scale differs from the correctly rounded prescribed rational.');
    end
    applied{k} = values.';
end
assert(isequal(sizes,double(entry.gram_sizes(:).')));
if isfield(raw,'gram_basis_exponents')
    bases = items(raw.gram_basis_exponents); assert(numel(bases) == 3);
    for k = 1:3
        E = double(bases{k}); assert(size(E,1) == sizes(k));
        derived = jet_targets(E,source,entry.q); exact = text_vector(targets{k});
        for j = 1:numel(derived), assert(qcompare(derived{j},qtext(exact{j})) == 0); end
    end
    basisStatus = 'Ordered exported Gram exponents checked against the exact manifold jet';
else
    assert(~mapping.historical_ordered_gram_basis_exponents_present);
    basisStatus = 'Not retained in the recovered historical files; no basis reconstructed';
end
coherence = struct('raw_copy',mapping.copy_path,'raw_sha256',mapping.sha256, ...
    'candidate_and_witness_match_summary',true,'Q_blocks_square_symmetric_and_finite',true, ...
    'gram_block_sizes',sizes,'actual_applied_scales',{applied}, ...
    'scaling_matches_variant_and_recorded_symbolic_targets',true, ...
    'ordered_gram_basis_status',basisStatus,'reported_reduced_gram_sizes',run.reduced_gram_sizes, ...
    'reduced_sizes_independently_reconstructed',false);
end

function [x,dW] = derivative(S,E,run)
n = double(S.N); x = qbinary(run.witness.x); ys = vector(run.witness.Y);
assert(numel(ys) == n); point = cell(n+1,1); point{1} = x;
for k = 1:n, point{k+1} = qbinary(ys(k)); end
a = flat_rationals(S.a); lam = flat_rationals(S.Lam); B = flat_rationals(S.B);
v = flat_rationals(S.v); w = flat_rationals(S.w); g = flat_rationals(S.g);
field = repmat({qinteger(0)},n+1,1);
for k = 1:numel(g), field{1} = qadd(field{1},qmul(g{k},qpower(x,k-1))); end
for i = 1:n
    field{1} = qadd(field{1},qmul(x,qmul(a{i,1},point{i+1})));
    for j = 1:n
        field{i+1} = qadd(field{i+1},qmul(lam{i,j},point{j+1}));
        field{i+1} = qadd(field{i+1},qmul(x,qmul(B{i,j},point{j+1})));
    end
    field{i+1} = qadd(field{i+1},qmul(qpower(x,2),v{i,1}));
    field{i+1} = qadd(field{i+1},qmul(qpower(x,3),w{i,1}));
end
E = double(E); coefficients = vector(run.candidate_coefficients);
assert(size(E,1) == numel(coefficients) && size(E,2) == n+1);
assert(all(E(:) >= 0 & E(:) == fix(E(:))));
gradient = repmat({qinteger(0)},n+1,1);
for i = 1:size(E,1)
    coefficient = qbinary(coefficients(i));
    for j = 1:n+1
        if E(i,j) == 0, continue; end
        term = qmul(coefficient,qinteger(E(i,j)));
        for k = 1:n+1
            term = qmul(term,qpower(point{k},E(i,k)-(k == j)));
        end
        gradient{j} = qadd(gradient{j},term);
    end
end
dW = qinteger(0);
for j = 1:n+1, dW = qadd(dW,qmul(gradient{j},field{j})); end
end

function targets = jet_targets(E,S,q)
assert(size(E,2) == S.N+1 && all(E(:) >= 0 & E(:) == fix(E(:))));
H = flat_rationals(S.H_exact); h = cell(S.N,1);
for j = 1:S.N
    h{j} = repmat({qinteger(0)},q+1,1);
    for degree = 1:q, h{j}{degree+1} = H{j,degree}; end
end
targets = cell(size(E,1),1);
for row = 1:size(E,1)
    term = repmat({qinteger(0)},q+1,1);
    if E(row,1) <= q, term{E(row,1)+1} = qinteger(1); end
    for j = 1:S.N
        for power = 1:E(row,j+1), term = truncated_product(term,h{j},q); end
    end
    top = qinteger(0);
    for degree = 0:q
        candidate = qabs(term{degree+1});
        if qcompare(candidate,top) > 0, top = candidate; end
    end
    if top.n.signum() == 0, targets{row} = qinteger(1);
    else, targets{row} = qdivide(qinteger(1),top); end
end
end

function product = truncated_product(a,b,q)
product = repmat({qinteger(0)},q+1,1);
for i = 0:q
    for j = 0:q-i
        product{i+j+1} = qadd(product{i+j+1},qmul(a{i+1},b{j+1}));
    end
end
end

function valid = rounds_to(target,value)
% Exact nearest-even rounding cell around the retained positive binary64.
assert(value > 0 && isfinite(value) && value < realmax);
bits = typecast(double(value),'uint64');
previous = typecast(bits-uint64(1),'double'); following = typecast(bits+uint64(1),'double');
middle = qbinary(value); two = qinteger(2);
lower = qdivide(qadd(qbinary(previous),middle),two);
upper = qdivide(qadd(middle,qbinary(following)),two);
left = qcompare(target,lower); right = qcompare(target,upper);
even = bitand(bits,uint64(1)) == 0;
valid = (left > 0 || (left == 0 && even)) && (right < 0 || (right == 0 && even));
end

function A = flat_rationals(field)
dims = double(field.dims(:)); values = text_vector(field.vals);
assert(numel(dims) == 2 && numel(values) == prod(dims));
A = cell(dims(1),dims(2));
for k = 1:numel(values), A{k} = qtext(values{k}); end
end

function q = qbinary(value)
assert(isscalar(value) && isreal(value) && isfinite(value));
bits = typecast(double(value),'uint64'); exponent = double(bitand(bitshift(bits,-52),uint64(2047)));
mantissa = bitand(bits,bitshift(uint64(1),52)-uint64(1));
if exponent == 0, shift = -1074;
else, mantissa = mantissa+bitshift(uint64(1),52); shift = exponent-1075; end
numerator = java.math.BigInteger(sprintf('%u',mantissa)); denominator = java.math.BigInteger('1');
if shift >= 0, numerator = numerator.shiftLeft(int32(shift));
else, denominator = denominator.shiftLeft(int32(-shift)); end
if bitget(bits,64), numerator = numerator.negate(); end
q = qmake(numerator,denominator);
end

function q = qtext(text)
text = strtrim(char(text)); pieces = strsplit(text,'/');
if numel(pieces) == 2
    q = qmake(java.math.BigInteger(pieces{1}),java.math.BigInteger(pieces{2})); return
end
assert(isscalar(pieces) && ~isempty(regexp(text, ...
    '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$','once')),'Invalid exact rational text.');
pieces = strsplit(lower(text),'e'); mantissa = pieces{1}; exponent = 0;
if numel(pieces) == 2, exponent = str2double(pieces{2}); end
dot = strfind(mantissa,'.'); decimals = 0;
if ~isempty(dot), decimals = numel(mantissa)-dot; end
mantissa = strrep(mantissa,'.',''); numerator = java.math.BigInteger(mantissa);
scale = exponent-decimals; ten = java.math.BigInteger('10'); denominator = java.math.BigInteger('1');
if scale >= 0, numerator = numerator.multiply(ten.pow(int32(scale)));
else, denominator = ten.pow(int32(-scale)); end
q = qmake(numerator,denominator);
end

function q = qinteger(value)
assert(isscalar(value) && isfinite(value) && value == fix(value));
q = struct('n',java.math.BigInteger(sprintf('%.0f',value)),'d',java.math.BigInteger('1'));
end

function q = qmake(n,d)
assert(d.signum() ~= 0,'Zero rational denominator.');
if d.signum() < 0, n = n.negate(); d = d.negate(); end
g = n.gcd(d); q = struct('n',n.divide(g),'d',d.divide(g));
end

function c = qadd(a,b)
if a.n.signum() == 0, c = b; return; elseif b.n.signum() == 0, c = a; return; end
c = qmake(a.n.multiply(b.d).add(b.n.multiply(a.d)),a.d.multiply(b.d));
end

function c = qmul(a,b)
if a.n.signum() == 0 || b.n.signum() == 0, c = qinteger(0); return; end
c = qmake(a.n.multiply(b.n),a.d.multiply(b.d));
end

function c = qdivide(a,b)
c = qmake(a.n.multiply(b.d),a.d.multiply(b.n));
end

function c = qpower(a,power)
assert(power >= 0 && power == fix(power));
c = struct('n',a.n.pow(int32(power)),'d',a.d.pow(int32(power)));
end

function c = qabs(a)
c = struct('n',a.n.abs(),'d',a.d);
end

function sign = qcompare(a,b)
sign = a.n.multiply(b.d).compareTo(b.n.multiply(a.d));
end

function text = qstring(a)
text = char(a.n.toString());
if ~a.d.equals(java.math.BigInteger('1')), text = [text '/' char(a.d.toString())]; end
end

function values = vector(value)
assert(isnumeric(value)); values = double(value(:));
end

function values = text_vector(value)
if ischar(value) || isstring(value), values = cellstr(value); return; end
assert(iscell(value)); values = cell(numel(value),1);
for k = 1:numel(value)
    current = value{k};
    while iscell(current) && isscalar(current), current = current{1}; end
    values{k} = char(current);
end
end

function rows = items(value)
if iscell(value), rows = value(:); else, rows = num2cell(value(:)); end
end

function check_hash(file,mapping)
info = dir(file); assert(isscalar(info) && info.bytes == mapping.bytes,'Retained file length changed.');
assert(strcmp(cm_sha256(file),mapping.sha256),'Retained file hash changed: %s',file);
end

function file = contained_path(root,relative)
base = char(java.io.File(root).getCanonicalPath());
file = char(java.io.File(fullfile(root,relative)).getCanonicalPath());
assert(startsWith(lower(file),[lower(base) filesep]),'Retained DD path escapes the reference root.');
end

function passed = arithmetic_self_checks()
% Synthetic checks do not infer any missing historical Gram basis.
one = qinteger(1); bits = typecast(1.0,'uint64'); next = typecast(bits+uint64(1),'double');
midpoint = qdivide(qadd(one,qbinary(next)),qinteger(2));
assert(rounds_to(midpoint,1.0) && ~rounds_to(midpoint,next),'Nearest-even tie check failed.');
assert(rounds_to(qtext('1/3'),1/3) && rounds_to(qtext('1/2'),0.5));
S = struct('N',1,'H_exact',struct('dims',[1,3],'vals',{{'0','2','3'}}));
targets = jet_targets([0,1;2,0;1,1;0,2],S,3);
expected = {'1/3','1','1/2','1'};
for k = 1:numel(targets), assert(qcompare(targets{k},qtext(expected{k})) == 0); end
passed = true;
end
