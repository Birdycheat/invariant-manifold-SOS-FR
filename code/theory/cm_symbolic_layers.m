function report = cm_symbolic_layers(selection)
% Exact restriction coefficients and moment functionals on every face direction.
report.curved_example = curved_example();
report.scalar_example = scalar_example();
cases=[2 3 2;2 5 3;3 3 2;3 5 3];
if strcmp(selection,'smoke'), cases=cases(2,:); end
records=cell(1,size(cases,1));
for k=1:size(cases,1), records{k}=layers(cases(k,1),cases(k,2),cases(k,3)); end
report.layers=records; report.passed=true;
fprintf('Exact exposing matrices: %d layer families and both curved examples\n',numel(records));
end

function out = layers(nx,m,degree)
E=up_to(nx+1,degree); alphas=up_to(nx,2*degree);
% R contains the coefficients of the Gram basis restricted to the graph.
R=sym(zeros(size(E,1),size(alphas,1)));
for i=1:size(E,1)
    powers=layer(nx,E(i,end));
    for j=1:size(powers,1)
        alpha=E(i,1:nx)+2*powers(j,:);
        [tf,id]=ismember(alpha,alphas,'rows'); assert(tf);
        R(i,id)=factorial(sym(E(i,end)))/prod(factorial(sym(powers(j,:))));
    end
end
q=floor(m/2); low=sum(alphas,2)<=q; C=R(:,low);
assert(rank(C)==nchoosek(nx+q,q));
records=cell(1,q+1);
for d=0:q
    % Rational box moments form the positive layer metric, not an identity matrix.
    ids=find(sum(alphas,2)==d); H=moments(alphas(ids,:),alphas(ids,:));
    for k=1:size(H,1), assert(isAlways(det(H(1:k,1:k))>0)); end
    Cd=R(:,ids); Z=Cd*H*Cd.'; assert(rank(Z)==rank(Cd));
    previous=R(:,sum(alphas,2)<d);
    if isempty(previous), B=sym(eye(size(E,1))); else, B=null(previous.'); end
    G=sym(zeros(size(E,1)));
    for r=0:2*d
        left=find(sum(alphas,2)==r); right=find(sum(alphas,2)==2*d-r);
        G=G+R(:,left)*moments(alphas(left,:),alphas(right,:))*R(:,right).';
    end
    % Test the exposing functional on every direction of the previous face.
    assert(all(isAlways(B.'*(G-Z)*B==0),'all'),'Affine layer functional differs');
    records{d+1}=struct('degree',d,'rank',double(rank(Z)), ...
        'C',{strings_of(Cd)},'H',{strings_of(H)},'Z',{strings_of(Z)}, ...
        'all_previous_face_directions_verified',true);
end
removed=E(:,end)==0 & sum(E(:,1:nx),2)<=q;
assert(nnz(removed)==rank(C));
out=struct('nx',nx,'m',m,'q',q,'basis_degree',degree,'gram_size',size(E,1), ...
    'eliminated_rank',double(rank(C)),'remaining_size',size(E,1)-double(rank(C)), ...
    'layers',{records},'adapted_deletion_exact',true,'passed',true);
end

function out = curved_example()
x1=sym('x1','real'); x2=sym('x2','real'); y=sym('y','real');
z=sym('z','real'); mu=sym('mu','positive');
phi=2*(x1^2+x2^2);
F=[-x1^5+x1*(y-phi);-x2^5+x2*(y-phi);-4*x1^6-4*x2^6+(2*phi-1)*(y-phi)];
assert(isequal(subs(jacobian(F,[x1;x2;y]),[x1;x2;y],[0;0;0]),diag(sym([0 0 -1]))));
moved=subs(F,y,z+phi); zdot=expand(moved(3)-4*x1*moved(1)-4*x2*moved(2));
assert(isAlways(zdot==-z));
C2=sym([0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1;2 0 2]);
alpha=[2 0;1 1;0 2]; H2=moments(alpha,alpha); Z=C2*H2*C2.'; Z2=Z(4:7,4:7);
assert(isequal(45*H2,sym([9 0 5;0 5 0;5 0 9])));
assert(isequal(45*Z2,sym([9 0 5 28;0 5 0 0;5 0 9 28;28 0 28 112])));
assert(all(isAlways(Z2*sym([-2;0;-2;1])==0)) && rank(Z2)==3);
for j=[1 3], M=C2(:,j)*C2(:,j).'; assert(any(diag(M)<sum(abs(M),2)-abs(diag(M)))); end
W=(x1^4+x2^4)/4+z^2/2;
dW=jacobian(W,[x1;x2;z])*[moved(1:2);zdot];
% Verify the explicit decrease certificate after the coordinate change.
target=(x1^4-z/2)^2+(x2^4-z/2)^2+z^2/2;
assert(isAlways(expand(-dW-target)==0));
Wq=(x1^2+x2^2)/2+mu*z^2;
dWq=jacobian(Wq,[x1;x2;z])*[moved(1:2);zdot];
assert(isAlways(simplify(-subs(dWq,[x2;z],[0;x1^2/(4*mu)])-x1^6+x1^4/(8*mu))==0));
out=struct('C2',{strings_of(C2)},'H2',{strings_of(H2)},'Z2',{strings_of(Z2)}, ...
    'kernel','y-2*x1^2-2*x2^2','adapted_field','(-x1^5+x1*z,-x2^5+x2*z,-z)', ...
    'certificate_identity',char(target),'quadratic_obstruction','x1^6-x1^4/(8*mu)','passed',true);
end

function out = scalar_example()
syms x y
phi=2*x^2+x^3; fx=x*(y-phi)-x^7; fy=-y+phi;
h=sym(0);
for d=2:10
    R=expand(subs(fy,y,h)-diff(h,x)*subs(fx,y,h)); c=coeffs(R,x,'All');
    if numel(c)>d, h=expand(h+c(end-d)*x^d); end
end
assert(isAlways(expand(h-2*x^2-x^3-4*x^8-3*x^9+16*x^10)==0));
E=up_to(2,3); C=sym(zeros(size(E,1),4));
for k=1:size(E,1)
    v=expand(x^E(k,1)*phi^E(k,2)); c=coeffs(v,x,'All');
    for d=0:min(3,numel(c)-1), C(k,d+1)=c(end-d); end
end
c2=sym(zeros(size(E,1),1)); c3=c2;
c2(ismember(E,[2 0],'rows'))=1; c2(ismember(E,[0 1],'rows'))=2;
c3(ismember(E,[3 0],'rows'))=1; c3(ismember(E,[0 1],'rows'))=1; c3(ismember(E,[1 1],'rows'))=2;
assert(isequal(C(:,3),c2) && isequal(C(:,4),c3));
assert(rank(C)==4);
for v={c2,c3}, M=v{1}*v{1}.'; assert(any(diag(M)<sum(abs(M),2)-abs(diag(M)))); end
Z=cell(1,4); for d=0:3, Z{d+1}=strings_of(C(:,d+1)*C(:,d+1).'/(2*d+1)); end
out=struct('jet_through_10',char(h),'basis_exponents',E,'C',{strings_of(C)}, ...
    'moment_exposers',{Z},'eliminated_rank',4,'remaining_size',size(E,1)-4, ...
    'mixed_vectors_not_dd',true,'passed',true);
end

function H = moments(A,B)
H=sym(zeros(size(A,1),size(B,1)));
for i=1:size(A,1)
    for j=1:size(B,1)
        power=A(i,:)+B(j,:);
        if all(mod(power,2)==0), H(i,j)=prod(1./sym(power+1)); end
    end
end
end

function E = up_to(n,d)
E=zeros(0,n); for k=0:d, E=[E;layer(n,k)]; end %#ok<AGROW>
end

function E = layer(n,d)
if n==1, E=d; return; end
E=zeros(0,n);
for j=d:-1:0, tail=layer(n-1,d-j); E=[E;j*ones(size(tail,1),1),tail]; end %#ok<AGROW>
end

function values = strings_of(P)
values=arrayfun(@char,P,'UniformOutput',false);
end
