function report = cm_symbolic_jets(root, selection)
% Reconstruct the finite invariant graph from the defining vector field.
benches = {'scalar','coupled','triple'};
folders = {fullfile(root,'data','benchmark'),fullfile(root,'data','coupled','instances'), ...
    fullfile(root,'data','triple','instances')};
records = {};
for b = 1:3
    files = dir(fullfile(folders{b},'*.json'));
    if strcmp(selection,'smoke')
        if b == 1
            names = {'N1_m5_s01.json','N3_m7_s01.json'};
        elseif b == 2
            names = {'cp_ny2_s01.json'};
        else
            names = {'tr_ny2_s01.json'};
        end
        files = files(ismember({files.name},names));
    end
    assert(~isempty(files),'No retained inputs found for %s',benches{b});
    for k = 1:numel(files)
        S = jsondecode(fileread(fullfile(files(k).folder,files(k).name)));
        [x,y,F,A,phi] = field_of(S,benches{b});
        N = S.m + 1;
        assert(all_zero(part(F(1:numel(x)),[x;y],1)), 'Center linear part is not zero');
        assert(all_zero(part(F(numel(x)+1:end)-A*y,[x;y],1)), 'Invalid stable linear part');
        assert(~isAlways(det(A)==0), 'Singular stable block');
        psi = sym(zeros(numel(y),1));
        snapshots = cell(1,N);
        % Cancel each homogeneous invariance residual using the stable block.
        for degree = 2:N
            restricted = expand(subs(F,y,psi));
            residual = restricted(numel(x)+1:end)-jacobian(psi,x)*restricted(1:numel(x));
            psi = expand(psi-A\homogeneous(residual,x,degree));
            snapshots{degree} = psi;
        end
        restricted = expand(subs(F,y,psi));
        residual = restricted(numel(x)+1:end)-jacobian(psi,x)*restricted(1:numel(x));
        red = part(restricted(1:numel(x)),x,S.m);
        checks = struct('invariance_residual_zero',all_zero(part(residual,x,N)), ...
            'reduced_order_exactly_m',all_zero(part(red,x,S.m-1)) && ~all_zero(red));
        if b == 1
            % Stored jets are comparison targets, not inputs to the recurrence.
            H = exact_flat(S.H_exact);
            expected = H(:,1:min(N,size(H,2)))*(x(1).^((1:min(N,size(H,2))).'));
            checks.stored_jet_matches = all_zero(psi-expected);
            checks.leading_field = all_zero(homogeneous(red,x,S.m)+x(1)^S.m);
        else
            eta = exact_matrix(S.eta4)*monomials(x,layer(numel(x),4));
            checks.T3_equals_field_shift = all_zero(snapshots{3}-phi);
            checks.T4_equals_shift_plus_eta = all_zero(snapshots{4}-phi-eta);
        end
        values = struct2cell(checks);
        rec = struct('bench',benches{b},'instance',files(k).name,'jet_order',N,'m',S.m, ...
            'checks',checks,'jet',{strings_of(psi)},'leading_field',{strings_of(homogeneous(red,x,S.m))}, ...
            'passed',all([values{:}]));
        records{end+1} = rec; %#ok<AGROW>
    end
    fprintf('Exact finite jets: %s, %d checked\n',benches{b},numel(files));
end
report = struct('count',numel(records),'records',{records},'passed',all(cellfun(@(r)r.passed,records)));
if strcmp(selection,'all'), assert(report.count==200,'Expected all 200 retained fields'); end
if strcmp(selection,'smoke'), assert(report.count==4,'Expected four smoke fields'); end
end

function [x,y,F,A,phi] = field_of(S,bench)
if strcmp(bench,'scalar')
    nx=1; ny=S.N; x=sym('x',[nx,1]); y=sym('y',[ny,1]);
    A=exact_flat(S.Lam); B=exact_flat(S.B); a=exact_flat(S.a);
    g=exact_flat(S.g); v=exact_flat(S.v); w=exact_flat(S.w);
    fx=g*(x(1).^((0:numel(g)-1).'))+x(1)*(a.'*y);
    F=[fx;A*y+x(1)*B*y+x(1)^2*v+x(1)^3*w]; phi=sym(zeros(ny,1));
    return
end
nx=S.nx; ny=S.ny; x=sym('x',[nx,1]); y=sym('y',[ny,1]);
A=exact_matrix(S.Lam); a=exact_matrix(S.a);
phi=exact_matrix(S.phi)*monomials(x,layer(nx,2)); z=y-phi;
f5=exact_matrix(S.Fhat)*monomials(x,layer(nx,5));
e4=exact_matrix(S.Ehat)*monomials(x,layer(nx,4));
stable=A;
for j=1:nx, stable=stable+x(j)*exact_matrix(S.(sprintf('B%d',j))); end
if strcmp(bench,'coupled')
    u=x(1)-x(2); rho=x.'*x; g1=sym(1); g2=sym(1);
    for j=1:floor((S.m+1)/2)
        g1=g1-(-1)^j*u^(2*j)/factorial(sym(2*j));
        g2=g2+(-rho)^j/factorial(sym(j));
    end
else
    g1=sym(1); g2=sym(1);
end
fx=g1*f5+x.*(a*z);
fz=stable*z+g2*e4;
F=part([fx;fz+jacobian(phi,x)*fx],[x;y],S.m+1);
end

function R = homogeneous(P,x,d)
t=sym('jet_degree_parameter'); R=sym(zeros(size(P)));
for k=1:numel(P)
    c=coeffs(expand(subs(P(k),x,t*x)),t,'All');
    if numel(c)>d, R(k)=c(end-d); end
end
end

function R = part(P,x,d)
t=sym('jet_degree_parameter'); R=sym(zeros(size(P)));
for k=1:numel(P)
    c=coeffs(expand(subs(P(k),x,t*x)),t,'All');
    last=min(d,numel(c)-1);
    for j=0:last, R(k)=R(k)+c(end-j); end
    R(k)=expand(R(k));
end
end

function E = layer(n,d)
if n==1, E=d; return; end
E=zeros(0,n);
for j=d:-1:0
    tail=layer(n-1,d-j); E=[E; j*ones(size(tail,1),1),tail]; %#ok<AGROW>
end
end

function v = monomials(x,E)
v=sym(ones(size(E,1),1));
for j=1:numel(x), v=v.*x(j).^E(:,j); end
end

function A = exact_matrix(value)
if iscell(value)
    if all(cellfun(@(q) ischar(q)||isstring(q),value(:)))
        A=sym(zeros(size(value))); for k=1:numel(value), A(k)=str2sym(value{k}); end
    else
        rows=cellfun(@(row) reshape(exact_matrix(row),1,[]),value(:),'UniformOutput',false); A=vertcat(rows{:});
    end
elseif isstring(value), A=str2sym(value);
elseif ischar(value), A=str2sym(value);
else, A=sym(value,'f');
end
end

function A = exact_flat(value)
A=exact_matrix(value.vals); A=reshape(A,double(value.dims(:).'));
end

function yes = all_zero(P)
yes=isequal(expand(P),sym(zeros(size(P))));
end

function values = strings_of(P)
values=arrayfun(@char,P,'UniformOutput',false);
end
