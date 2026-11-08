function [report, model] = cm_beg_bounds(S)
% Reconstruct T5/T11 and prove every stored BEG Cauchy enclosure independently.
syms x y z
    a=q(S.a); b=-a*(a-3)/(2*(a-1)); lam=3*(a-1)/(2*a); h2=a/3;
    c5=a^4*(a^2-2*a-9)/(270*(a-1)^2);
    assert(isAlways(q(S.b)==b & q(S.rate)==lam & q(S.h2)==h2 & q(S.c5)==c5));
    F5=field_xy(5,a,b,x,y); h=sym(0);
    % Solve the invariance equation degree by degree from the field coefficients.
    for d=2:9
        residual=expand(subs(F5(2),y,h)-diff(h,x)*subs(F5(1),y,h));
        h=expand(h+coefficient(residual,x,d)*x^d/lam);
    end
    reduced=part(subs(F5(1),y,h),x,9);
    assert(isAlways(coefficient(h,x,2)==h2));
    assert(isAlways(part(reduced,x,4)==0 & coefficient(reduced,x,5)==c5));
    assert_zero(h-series_from_struct(S.h_xy,x));
    assert_zero(reduced-series_from_struct(S.reduced,x));
    SX=expand(subs(F5,y,z+h2*x^2)); SX(2)=expand(SX(2)-2*h2*x*SX(1));
    F11=field_xy(11,a,b,x,y);
    T=expand(subs(F11,y,z+h2*x^2)); T(2)=expand(T(2)-2*h2*x*T(1));
    T=part(T,[x;z],11); T5=part(T,[x;z],5);
    assert_zero(F5-fields_from_rows(S.Sxy,x,y));
    assert_zero(SX-fields_from_rows(S.Sxz,x,z));
    assert_zero(T5-fields_from_rows(S.T5xz,x,z));
    assert_zero(T-fields_from_rows(S.T11xz,x,z));
    coefficientNorms=sym(zeros(2,6));
    for i=1:2
        for d=6:11, coefficientNorms(i,d-5)=sum(abs(coeffs(homogeneous(T(i),[x;z],d),[x;z]))); end
    end
    % Search admissible complex polydiscs for a tighter rational Cauchy bound.
    grid=cauchy_grid(a,b); bounds=cell(1,numel(S.radii));
    exactBounds=struct('r',{},'delta',{},'S',{},'B',{},'M_i',{},'M',{});
    for j=1:numel(S.radii)
        r=q(S.radii(j).r); [B,delta,bigS,id]=best_bound(grid,r);
        % Combine explicit degrees 6--11 with the analytic degree-12 tail.
        Mi=coefficientNorms*(r.^((0:5).'))+B*r^6;
        M=Mi(1); if isAlways(Mi(2)>M), M=Mi(2); end
        storedMi=q_vector(S.radii(j).M_i);
        assert_zero(Mi-storedMi);
        assert(isAlways(B==q(S.radii(j).B) & delta==q(S.radii(j).delta) & M==q(S.radii(j).M)));
        u=delta+a*delta^2/3; v=abs(b)*u; eta=grid.ex(id)*grid.ch(id)-1;
        gap=1+1/(a-1)-eta/(a-1);
        assert(isAlways(v<16 & (a*delta)^2<23*24 & (a*delta)^2<24*25));
        assert(isAlways(gap>0 & eta<a & 1-sym(3)*r/(2*delta)>0));
        assert(isAlways(B==64*bigS/delta^12/(1-sym(3)*r/(2*delta))));
        assert(isAlways(M>0 & M>=Mi(1) & M>=Mi(2)));
        % Keep the certified real ball inside the physical spin domain.
        assert(isAlways(1/a-r>r & 1/a+r+a*r^2/3<1));
        exactBounds(j)=struct('r',r,'delta',delta,'S',bigS,'B',B,'M_i',Mi,'M',M);
        bounds{j}=struct('r',char(r),'complex_polydisc_radius',char(delta), ...
            'denominator_lower_bound',char(gap),'S',char(bigS),'B',char(B), ...
            'M_i',{strings_of(Mi)},'M',char(M),'exact_reference_match',true,'physical_domain_verified',true);
        
    end
    if isAlways(c5<0), classification='stable'; else, classification='unstable'; end
    assert(strcmp(classification,S.class));
    report=struct('status','passed','instances',1,'bounds_checked',numel(bounds), ...
        'T5_checked',true,'T11_checked',true,'a',char(a),'class',classification, ...
        'jet_through_9',char(h),'reduced_surrogate_through_9',char(reduced), ...
        'field_coefficients_reconstructed',true,'bounds',{bounds},'passed',true);
model=struct('a',a,'b',b,'rate',lam,'h2',h2,'x',x,'z',z, ...
    'T5xz',T5,'T11xz',T,'bounds',exactBounds);
end

function F = field_xy(N,a,b,x,y)
rho=1/(a-1); C=sym(zeros(1,N+1)); A=C;
for d=0:N
    for j=0:d
        k=d-j; coefficient=rho*b^j*a^k/(factorial(sym(j))*factorial(sym(k)));
        if mod(k,2)==0, C(d+1)=C(d+1)+coefficient*x^k*y^j;
        else, A(d+1)=A(d+1)+coefficient*x^k*y^j; end
    end
end
G=sym(zeros(2,N+1)); denominator=1+rho;
for d=0:N
    term=[A(d+1);C(d+1)];
    for j=1:d, term=term-C(j+1)*G(:,d-j+1); end
    G(:,d+1)=expand(term/denominator);
end
F=expand(sum(G,2)+[-x;-1/a-y]);
assert_zero(subs(F,[x;y],[0;0]));
end

function grid = cauchy_grid(a,b)
d=sym((4:40).')/100; u=d+a*d.^2/3; v=abs(b)*u; ad=a*d;
% Rational Taylor sums plus geometric tails bound exp, cosh and sinh.
ex=sym(zeros(size(d))); ch=ex; sh=ex;
for k=0:14, ex=ex+v.^k/factorial(sym(k)); end
ex=ex+v.^15/factorial(sym(15))./(1-v/16);
for k=0:10
    ch=ch+ad.^(2*k)/factorial(sym(2*k));
    sh=sh+ad.^(2*k+1)/factorial(sym(2*k+1));
end
ch=ch+ad.^22/factorial(sym(22))./(1-ad.^2/(23*24));
sh=sh+ad.^23/factorial(sym(23))./(1-ad.^2/(24*25));
gap=1+1/(a-1)-(ex.*ch-1)/(a-1);
sx=d+ex.*sh./((a-1)*gap); sy=1/a+u+ex.*ch./((a-1)*gap);
sz=sy+2*a*d.*sx/3; S=sx;
for k=1:numel(d), if isAlways(sz(k)>S(k)), S(k)=sz(k); end, end
grid=struct('d',d,'ex',ex,'ch',ch,'gap',gap,'S',S);
end

function [best,d,S,index] = best_bound(grid,r)
index=0; best=sym(0); d=sym(0); S=sym(0);
for k=1:numel(grid.d)
    delta=grid.d(k);
    if ~isAlways(sym(3)*r/(2*delta)<sym(1)/2 && grid.gap(k)>0), continue; end
    value=64*grid.S(k)/delta^12/(1-sym(3)*r/(2*delta));
    if index==0 || isAlways(value<best)
        best=value; d=delta; S=grid.S(k); index=k;
    end
end
assert(index>0,'No admissible Cauchy polydisc');
end

function F = fields_from_rows(rows,x,w)
F=sym(zeros(2,1));
if iscell(rows) && size(rows,2)==5, count=size(rows,1); matrixRows=true;
else, count=numel(rows); matrixRows=false; end
for k=1:count
    if matrixRows, row=rows(k,:); else, row=rows{k}; end
    component=double(row{1}); i=double(row{2}); j=double(row{3});
    F(component)=F(component)+q(row{4})*x^i*w^j;
end
F=expand(F);
end

function P = series_from_struct(S,x)
P=sym(0); names=fieldnames(S);
for k=1:numel(names)
    exponent=str2double(regexp(names{k},'\d+','match','once'));
    P=P+q(S.(names{k}))*x^exponent;
end
end

function R = homogeneous(P,x,d)
t=sym('beg_degree_parameter'); R=sym(zeros(size(P)));
for k=1:numel(P), c=coeffs(expand(subs(P(k),x,t*x)),t,'All'); if numel(c)>d, R(k)=c(end-d); end, end
end

function R = part(P,x,d)
R=sym(zeros(size(P))); for j=0:d, R=R+homogeneous(P,x,j); end
R=expand(R);
end

function c = coefficient(P,x,d)
all=coeffs(expand(P),x,'All'); c=sym(0); if numel(all)>d, c=all(end-d); end
end

function value = q(value)
if iscell(value), value=value{1}; end
if ischar(value)||isstring(value), value=str2sym(value); else, value=sym(value,'f'); end
end

function v = q_vector(values)
v=sym(zeros(numel(values),1)); for k=1:numel(values), if iscell(values), item=values{k}; else, item=values(k); end, v(k)=q(item); end
end

function assert_zero(P)
assert(all(isAlways(simplify(expand(P))==0),'all'),'Exact polynomial comparison failed');
end

function values = strings_of(P)
values=arrayfun(@char,P,'UniformOutput',false);
end
