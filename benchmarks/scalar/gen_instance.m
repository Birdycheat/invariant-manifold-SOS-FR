function P = gen_instance(N, m, d, seed, K, wscale, ascale)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if nargin < 5 || isempty(K), K = 24; end
if nargin < 6 || isempty(wscale), wscale = 1; end
if nargin < 7, ascale = 1; end
rng(seed);
lam = 1 + 2*rand(N,1);
Lam = -diag(lam);
B = 0.8*(2*rand(N) - 1);
a = randn(N,1); a = ascale*a/norm(a);
v = lam.*a + 0.3*randn(N,1);
v = v/(a'*(-Lam\v));
w = wscale*randn(N,1);
g = zeros(1, m+1); g(4) = -1;
H = zeros(N, K);
H(:,2) = -Lam\v;
% Cancel lower reduced terms and retain the prescribed -d*x^m central term.
for k = 3:K
    if k >= 4 && k <= m
        target = -d*(k == m);
        g(k+1) = target - a'*H(:,k-1);
    end
    acc = -B*H(:,k-1) - w*(k == 3);
    for j = 2:k-2
        acc = acc + j*H(:,j)*Fk(k+1-j, g, a, H, m);
    end
    H(:,k) = Lam\acc;
end
P = struct('N',N,'m',m,'d',d,'seed',seed,'lam',lam,'Lam',Lam,'B',B,'a',a,'v',v,'w',w,'g',g,'H',H);
end

function c = Fk(k, g, a, H, m)
c = a'*H(:,k-1);
if k <= m, c = c + g(k+1); end
end
