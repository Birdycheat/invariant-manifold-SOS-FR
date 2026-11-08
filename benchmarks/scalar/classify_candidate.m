function [st, rb, cand] = classify_candidate(ok, info, P, r, p)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
st = 'echec solveur'; rb = NaN; cand = struct('E', [], 'c', [], 'x', NaN, 'Y', NaN);
if ~ok, return; end
E = info.E; c = info.cval;
fx = @(X, Y) polyval(fliplr(P.g), X) + X.*(P.a'*Y);
fy = @(X, Y) P.Lam*Y + (P.B*Y).*X + P.v*X.^2 + P.w*X.^3;
if p > 0
    Hp = P.H(:,2:p);
    phiX  = @(X) Hp*(X.^((2:p)'));
    dphiX = @(X) Hp*((2:p)'.*X.^((1:p-1)'));
    dWn = @(X, Y) sum(gradpolyE(E, c, [X; Y - phiX(X)]).*[fx(X, Y); fy(X, Y) - dphiX(X).*fx(X, Y)], 1);
else
    dWn = @(X, Y) sum(gradpolyE(E, c, [X; Y]).*[fx(X, Y); fy(X, Y)], 1);
end
dWe = @(x, Y) exact_dW(E, c, P, x, Y, p);
% Use numerical searches to propose points, then check their derivative signs exactly.
[st, rb, xb, yb] = exact_classify(dWn, dWe, P, r);
cand = struct('E', E, 'c', c, 'x', xb, 'Y', yb);
end

function [status, rb, xbad, ybad] = exact_classify(dWnum, dWex, P, r)
[rb, ~, pt] = validate_nearN(dWnum, P, r);
xbad = NaN; ybad = NaN(P.N, 1);
if rb > 0 && isAlways(dWex(pt.x, pt.Y) > 0)
    status = 'invalide'; xbad = pt.x; ybad = pt.Y; return
end
rb = 0;
m = P.m; N = P.N;
xmin = floor_x(m);
xg = logspace(log10(xmin), -9, round(log10(xmin)) + 10);
xg = reshape([xg; -xg], 1, []);
K = size(P.H, 2); hx = @(X) P.H(:,2:K)*(X.^((2:K)'));
if N == 1, dirs = 1; else, dirs = [eye(N), P.a/norm(P.a)]; end
for xx = xg
    for q = 2:4
        for t = [-1 -0.1 0.1 1]
            for e = dirs
                Y = hx(xx) + t*abs(xx)^q*e;
                if isAlways(dWex(xx, Y) > 0)
                    status = 'invalide'; xbad = xx; ybad = Y; rb = abs(xx); return
                end
            end
        end
    end
end
status = 'aucune violation';
end

function [rbad, mx, pt] = validate_nearN(dWfun, P, r)
N = P.N; m = P.m;
xmin = floor_x(m);
xs = [logspace(log10(xmin), log10(r), 80), -logspace(log10(xmin), log10(r), 80)];
rs = RandStream('mt19937ar', 'Seed', 7);
Dirs = [eye(N), -eye(N), P.a, -P.a, randn(rs, N, 8)];
Dirs = Dirs./vecnorm(Dirs);
t = linspace(-3, 3, 25);
K = size(P.H, 2);
hx = @(X) P.H(:,2:K)*(X.^((2:K)'));
colmax = -inf(numel(xs), 1); mx = -inf;
[XI, TI, DI] = ndgrid(1:numel(xs), 1:numel(t), 1:size(Dirs,2));
X = xs(XI(:));
colpt = cell(numel(xs), 1);
for q = 1:4
    Z = Dirs(:, DI(:)).*(t(TI(:)).*abs(X).^q);
    Y = hx(X) + Z;
    Dv = dWfun(X, Y);
    mx = max(mx, max(Dv./(abs(X).^(2*m-2) + sum(Z.^2, 1))));
    for ix = 1:numel(xs)
        k = find(XI(:) == ix);
        [v, j] = max(Dv(k));
        if v > colmax(ix)
            colmax(ix) = v;
            colpt{ix} = struct('x', X(k(j)), 'Y', Y(:, k(j)), 'val', v);
        end
    end
end
bad = find(colmax > 0);
pt = [];
if isempty(bad)
    rbad = 0;
else
    [rbad, i] = max(abs(xs(bad)));
    pt = colpt{bad(i)};
end
end

function xmin = floor_x(m)
if m > 3, xmin = 10^(-12/(m-3)); else, xmin = 1e-9; end
end

function v = exact_dW(E, c, P, x, Y, p)
% Preserve the binary64 values defining the floating candidate and coordinate shift.
S = @(a) sym(a, 'f');
[Lam, B, a, vv, w] = exact_data(P);
xs = S(x); Ys = S(Y(:));
g0 = sum(P.gS.*xs.^(0:numel(P.gS)-1));
fx = g0 + xs*(a.'*Ys);
fy = Lam*Ys + xs*(B*Ys) + vv*xs^2 + w*xs^3;
if p > 0
    Hp = S(P.H(:,2:p));
    phi = Hp*(xs.^((2:p).')); dphi = Hp*((2:p).'.*xs.^((1:p-1).'));
    pts = [xs; Ys - phi];
    fv = [fx; fy - dphi*fx];
else
    pts = [xs; Ys];
    fv = [fx; fy];
end
keep = c ~= 0;
E = E(keep,:); cs = S(c(keep)); cs = cs(:);
nv = size(E, 2);
G = sym(zeros(nv, 1));
for j = 1:nv
    k = E(:,j) > 0;
    if ~any(k), continue; end
    Ej = E(k,:); Ej(:,j) = Ej(:,j) - 1;
    mon = prod(repmat(pts.', size(Ej,1), 1).^Ej, 2);
    G(j) = sum(cs(k).*E(k,j).*mon);
end
v = G.'*fv;
end

function G = gradpolyE(E, c, Pts)
[nm, nv] = size(E); K = size(Pts, 2);
G = zeros(nv, K);
for j = 1:nv
    acc = zeros(1, K);
    for i = 1:nm
        if E(i,j) == 0 || c(i) == 0, continue; end
        term = c(i)*E(i,j)*ones(1, K);
        for k = 1:nv
            e = E(i,k) - (k == j);
            if e > 0, term = term.*Pts(k,:).^e; end
        end
        acc = acc + term;
    end
    G(j,:) = acc;
end
end
