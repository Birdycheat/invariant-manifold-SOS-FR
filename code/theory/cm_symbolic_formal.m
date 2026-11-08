function report = cm_symbolic_formal()
% Exact algebraic examples delimiting the facial-reduction statements.
x=sym('x','real'); y=sym('y','real'); x1=sym('x1','real'); x2=sym('x2','real');
checks={};

W=(x^2+y^2)/2; F=[-x^5;-y]; v=[x^3;y]; Q=sym(eye(2))/2;
zero(-jacobian(W,[x;y])*F-(x^6+y^2)/2-v.'*Q*v);
zero(subs(jacobian(F,[x;y]),[x;y],[0;0])-diag(sym([0 -1])));
checks{end+1}=struct('name','nonhyperbolic_strict_gram','Q',{text(Q)},'restricted_jet_through_2','0','passed',true);

W=x^2/2; gamma=x^2; sigma=x^2/2;
zero(-diff(W,x)*(-x^3)-x^4/2-gamma*sigma);
assert(isAlways(subs(gamma,x,0)==0));
checks{end+1}=struct('name','vanishing_localizer','gamma',char(gamma),'sigma',char(sigma), ...
    'positive_gram','1/2','degree_one_kernel_forcing',false,'passed',true);

% The moment functional must vanish on every Gram representation of zero.
v=[x1^2;x1*x2;x2^2]; B=sym([0 0 1;0 -2 0;1 0 0]);
H=[sym(1)/5 0 sym(1)/9;0 sym(1)/9 0;sym(1)/9 0 sym(1)/5];
zero(v.'*B*v); zero(trace(H*B)); assert(trace(B)==-2);
assert(isAlways(det(H)==sym(56)/18225));
for k=1:3, assert(isAlways(det(H(1:k,1:k))>0)); end
assert(H(1,3)==H(2,2));
checks{end+1}=struct('name','moments_preserve_polynomial_identities','H',{text(H)}, ...
    'null_polynomial_gram',{text(B)},'moment_pairing','0','identity_pairing','-2','passed',true);

v=[sym(1);x;x^2;y]; Q=diag(sym([0 -1 2 2])); Z=sym(zeros(4)); Z(2,2)=sym(1)/3;
zero(-jacobian(x^2+y^2,[x;y])*[-x^3;-y]-x^2-v.'*Q*v);
zero(Q(:,1)); assert(isAlways(trace(Z*Q)==-sym(1)/3));
checks{end+1}=struct('name','margin_order_is_required','Q_affine',{text(Q)}, ...
    'degree_one_exposer_pairing','-1/3','passed',true);

% Diagonal basis scaling preserves the polynomial but can change DD membership.
c=sym([1;2]); D=diag([sym(1),sym(1)/2]); v=[x^2;y]; cs=D*c;
aa=sym('aa','real'); bb=sym('bb','real'); dd=sym('dd','real');
Q=[aa bb;bb dd]; scaled=D.'\Q/D;
zero(v.'*Q*v-(D*v).'*scaled*(D*v)); zero(det(scaled)-4*det(Q));
zero(trace((c*c.')*Q)-trace((cs*cs.')*scaled));
Q0=sym([4 -2;-2 1]); scaled0=D.'\Q0/D;
zero(Q0*c); zero(scaled0*cs); zero(v.'*Q0*v-(2*x^2-y)^2);
assert(~dd_matrix(c*c.') && dd_matrix(cs*cs.') && ~dd_matrix(Q0) && dd_matrix(scaled0));
checks{end+1}=struct('name','dd_depends_on_diagonal_basis_scaling','Q',{text(Q0)}, ...
    'scaled_Q',{text(scaled0)},'same_polynomial','(2*x^2-y)^2','passed',true);

v=[x;x^2]; decay=x^4;
zero(decay-v.'*diag(sym([0 1]))*v); zero(decay+x^2-v.'*sym(eye(2))*v);
zero(-diff(-x^2/2,x)*x^3-x^4/2-(x^2)^2/2);
checks{end+1}=struct('name','full_feasible_set_includes_candidate_and_margin_constraints', ...
    'excluded_margins',{{'0','-1'}},'decrease_sos_does_not_impose_W_positive',true,'passed',true);

% Rational field coefficients alone do not ensure a selected graph is rational.
parameter=sym('parameter','real'); F=[-x^5;-2*x^4*y]; h=parameter*x^2;
zero(subs(F(2),y,h)-diff(h,x)*F(1)); zero(subs(jacobian(F,[x;y]),[x;y],[0;0]));
c=[sym(1);sqrt(sym(2))]; Z=c*c.'/5; t=sym('t');
zero(5^2*Z(1,2)^2-2); zero(det(Z)); assert(rank(Z)==1);
checks{end+1}=struct('name','rational_field_can_have_an_irrational_invariant_graph', ...
    'field',{{char(F(1)),char(F(2))}},'graph','sqrt(2)*x^2','C',{text(c)}, ...
    'computed_exposer',{text(Z)},'off_diagonal_minimal_polynomial',char(25*t^2-2), ...
    'invertible_transverse_linear_block',false,'passed',true);

C=sym([1 0 0;0 1 0;0 0 1;2 0 2]); Z=sym(3)*C*H*C.'/2;
kernel=sym([-2;0;-2;1]); zero(Z*kernel); assert(rank(Z)==3);
checks{end+1}=struct('name','rational_jet_and_localizer_give_rational_exposers', ...
    'gamma_at_origin','3/2','Z',{text(Z)},'rational_kernel',{text(kernel)},'passed',true);
report=struct('checks_passed',numel(checks),'checks',{checks},'passed',true);
fprintf('Exact formal edge cases: %d checked\n',numel(checks));
end

function zero(P)
assert(all(isAlways(simplify(expand(P))==0),'all'),'Exact identity failed');
end

function yes = dd_matrix(M)
yes=all(isAlways(diag(M)>=sum(abs(M),2)-abs(diag(M))));
end

function values = text(A)
values=arrayfun(@char,A,'UniformOutput',false);
end
