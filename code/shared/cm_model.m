function [fields, phi, meta] = cm_model(bench, data, coordinates, option)
%CM_MODEL Reconstruct exact polynomial fields from the instance coefficients.
if nargin < 3, coordinates = 'xz'; end
if nargin < 4, option = []; end
switch lower(char(bench))
    case 'scalar'
        [fields, phi, meta] = scalarModel(data, coordinates, option);
    case 'coupled'
        if isempty(option), option = 1; end
        [fields, phi, meta] = coupledModel(data, coordinates, option);
    case 'triple'
        [fields, phi, meta] = tripleModel(data, coordinates);
    case 'beg'
        [fields, phi, meta] = begModel(data, coordinates, option);
    otherwise
        error('cm:model:bench', 'Unknown benchmark: %s', bench);
end
end

function [fields, phi, meta] = scalarModel(S, coordinates, option)
ny = double(S.N); nv = ny+1; X = variables(nv);
Lam = cm_exact_flat(S.Lam); B = cm_exact_flat(S.B);
a = cm_exact_flat(S.a); v = cm_exact_flat(S.v); w = cm_exact_flat(S.w);
g = cm_exact_flat(S.g); g = g(:);
if isstruct(option) && isfield(option,'H')
    H = cm_exact(option.H);
elseif isfield(S,'H')
    H = cm_exact_flat(S.H);
else
    H = cm_exact_flat(S.H_exact);
end
if isstruct(option) && isfield(option,'p'), H = H(:,1:double(option.p)); end
phi = cell(1,ny);
for k = 1:ny
    phi{k} = cm_poly('zero',nv);
    for j = 1:size(H,2), phi{k} = add(phi{k}, monomial(nv,1,j), H(k,j)); end
end
% Substitute y=z+phi; the transverse equation also needs the chain rule.
y = X(2:end);
if strcmp(coordinates,'xz')
    for k = 1:ny, y{k} = add(y{k},phi{k}); end
else
    assert(strcmp(coordinates,'xy'), 'cm:model:coordinates', 'Expected xy or xz.');
end
fx = cm_poly('zero',nv);
for j = 1:numel(g), fx = add(fx,monomial(nv,1,j-1),g(j)); end
for k = 1:ny, fx = add(fx,mul(X{1},y{k}),a(k)); end
f = [{fx},cell(1,ny)];
for k = 1:ny
    fy = add(cm_poly('scale',monomial(nv,1,2),v(k)),monomial(nv,1,3),w(k));
    for j = 1:ny
        fy = add(fy,y{j},Lam(k,j));
        fy = add(fy,mul(X{1},y{j}),B(k,j));
    end
    if strcmp(coordinates,'xz'), fy = add(fy,mul(cm_poly('deriv',phi{k},1),fx),-1); end
    f{k+1} = fy;
end
fields = {f}; meta = struct('nx',1,'ny',ny,'coordinates',coordinates);
end

function [fields, phi, meta] = coupledModel(S, coordinates, enclosure)
ny = double(S.ny); nv = ny+2; X = variables(nv);
a = cm_exact(S.a); Lam = cm_exact(S.Lam); B = {cm_exact(S.B1),cm_exact(S.B2)};
F = cm_exact(S.Fhat); E = cm_exact(S.Ehat); P = cm_exact(S.phi); H = cm_exact(S.eta4);
phi = forms(P,2,nv); Fh = forms(F,5,nv); Eh = forms(E,4,nv); eta = forms(H,4,nv);
z = X(3:end);
for k = 1:ny
    if strcmp(coordinates,'xy'), z{k} = add(z{k},phi{k},-1);
    elseif strcmp(coordinates,'xz4'), z{k} = add(z{k},eta{k});
    else, assert(strcmp(coordinates,'xz'), 'cm:model:coordinates', 'Unknown coordinates.');
    end
end
f0 = cell(1,nv); f1 = cell(1,nv); f2 = cell(1,nv);
for i = 1:2
    f0{i} = cm_poly('zero',nv); f1{i} = Fh{i}; f2{i} = cm_poly('zero',nv);
    for k = 1:ny, f0{i} = add(f0{i},mul(X{i},z{k}),a(i,k)); end
end
for k = 1:ny
    f0{k+2} = cm_poly('zero',nv); f1{k+2} = cm_poly('zero',nv); f2{k+2} = Eh{k};
    for j = 1:ny
        f0{k+2} = add(f0{k+2},z{j},Lam(k,j));
        for i = 1:2, f0{k+2} = add(f0{k+2},mul(X{i},z{j}),B{i}(k,j)); end
    end
    if ~strcmp(coordinates,'xz')
        if strcmp(coordinates,'xy'), shape = phi{k}; direction = 1;
        else, shape = eta{k}; direction = -1; end
        for i = 1:2
            derivative = cm_poly('deriv',shape,i);
            f0{k+2} = add(f0{k+2},mul(derivative,f0{i}),direction);
            f1{k+2} = add(f1{k+2},mul(derivative,Fh{i}),direction);
        end
    end
end
one = cm_poly('mono',zeros(1,nv)); u = add(X{1},X{2},-1);
rho = add(mul(X{1},X{1}),mul(X{2},X{2})); u2 = mul(u,u);
g1hi = add(one,u2,sym(1)/2); g2lo = add(one,rho,-1);
if enclosure == 1
    G1 = {one,g1hi}; G2 = {g2lo,one};
elseif enclosure == 2
    G1 = {add(g1hi,mul(u2,u2),-sym(1)/24),g1hi};
    G2 = {g2lo,add(g2lo,mul(rho,rho),sym(1)/2)};
else
    error('cm:model:enclosure', 'Only the proved enclosure orders 1 and 2 are supported.');
end
% The two independent analytic intervals give four polynomial vertices.
fields = cell(1,4); vertex = 0;
for i = 1:2
    for j = 1:2
        vertex = vertex+1; f = cell(1,nv);
        for k = 1:nv, f{k} = add(add(f0{k},mul(G1{i},f1{k})),mul(G2{j},f2{k})); end
        fields{vertex} = f;
    end
end
meta = struct('nx',2,'ny',ny,'coordinates',coordinates, ...
    'enclosure_order',enclosure,'analytic_enclosures','global alternating Taylor inequalities');
meta.parts = {f0,f1,f2}; meta.enclosures = {G1,G2};
end

function [fields, phi, meta] = tripleModel(S, coordinates)
nx = 3; ny = double(S.ny); nv = nx+ny; X = variables(nv);
a = cm_exact(S.a); Lam = cm_exact(S.Lam);
B = {cm_exact(S.B1),cm_exact(S.B2),cm_exact(S.B3)};
phi = forms3(cm_exact(S.phi),2,nv); F = forms3(cm_exact(S.Fhat),5,nv); E = forms3(cm_exact(S.Ehat),4,nv);
z = X(nx+1:end);
if strcmp(coordinates,'xy')
    for k = 1:ny, z{k} = add(z{k},phi{k},-1); end
else
    assert(strcmp(coordinates,'xz'), 'cm:model:coordinates', 'Expected xy or xz.');
end
f = cell(1,nv);
for i = 1:nx
    f{i} = F{i};
    for k = 1:ny, f{i} = add(f{i},mul(X{i},z{k}),a(i,k)); end
end
for k = 1:ny
    f{nx+k} = E{k};
    for j = 1:ny
        f{nx+k} = add(f{nx+k},z{j},Lam(k,j));
        for i = 1:nx, f{nx+k} = add(f{nx+k},mul(X{i},z{j}),B{i}(k,j)); end
    end
    if strcmp(coordinates,'xy')
        for i = 1:nx, f{nx+k} = add(f{nx+k},mul(cm_poly('deriv',phi{k},i),f{i})); end
    end
end
fields = {f}; meta = struct('nx',nx,'ny',ny,'coordinates',coordinates);
end

function [fields, phi, meta] = begModel(S, coordinates, option)
nv = 2; phi = {cm_poly('mono',[2,0],cm_exact(S.h2))};
if isstruct(option)
    assert(strcmp(coordinates,'xz'));
    rows = cm_json_rows(S.T5xz); base = rowFields(rows,nv);
    rho = add(monomial(nv,1,2),monomial(nv,2,2));
    % The certified BEG remainder is bounded componentwise by M*(x^2+z^2)^3.
    bound = cm_poly('scale',mul(mul(rho,rho),rho),cm_exact(option.M));
    fields = cell(1,4); vertex = 0;
    for a = [-1,1]
        for b = [-1,1]
            vertex = vertex+1; fields{vertex} = {add(base{1},bound,a),add(base{2},bound,b)};
        end
    end
elseif strcmp(coordinates,'xy')
    fields = {rowFields(cm_json_rows(S.Sxy),nv)};
else
    assert(strcmp(coordinates,'xz'));
    fields = {rowFields(cm_json_rows(S.Sxz),nv)};
end
meta = struct('nx',1,'ny',1,'coordinates',coordinates);
end

function out = rowFields(rows,nv)
out = {cm_poly('zero',nv),cm_poly('zero',nv)};
for k = 1:size(rows,1)
    component = double(cm_exact(rows{k,1}));
    exponent = [double(cm_exact(rows{k,2})),double(cm_exact(rows{k,3}))];
    out{component} = add(out{component},cm_poly('mono',exponent,cm_exact(rows{k,4})));
end
end

function out = forms(matrix,degree,nv)
out = cell(1,size(matrix,1)); exponents = [(degree:-1:0).',(0:degree).',zeros(degree+1,nv-2)];
for k = 1:numel(out), out{k} = cm_poly('from',exponents,matrix(k,:).'); end
end

function out = forms3(matrix,degree,nv)
exponents = zeros((degree+1)*(degree+2)/2,nv); position = 0;
for a = degree:-1:0
    for b = degree-a:-1:0, position = position+1; exponents(position,1:3) = [a,b,degree-a-b]; end
end
out = cell(1,size(matrix,1));
for k = 1:numel(out), out{k} = cm_poly('from',exponents,matrix(k,:).'); end
end

function out = variables(nv)
out = cell(1,nv);
for j = 1:nv, out{j} = monomial(nv,j,1); end
end

function value = monomial(nv,j,power)
exponent = zeros(1,nv); exponent(j) = power; value = cm_poly('mono',exponent);
end

function value = add(a,b,varargin)
value = cm_poly('add',a,b,varargin{:});
end

function value = mul(a,b)
value = cm_poly('mul',a,b);
end
