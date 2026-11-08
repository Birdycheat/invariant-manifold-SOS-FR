function [gS, HS] = exact_jet(P, Kmax, tune)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
if nargin < 3, tune = true; end
N = P.N; m = P.m;
[Lam, B, a, v, w] = exact_data(P);
if tune
    gS = sym(zeros(1, m+1));
else
    gS = sym(P.g, 'f');
end
HS = sym(zeros(N, Kmax));
HS(:,2) = -Lam\v;
% Solve the invariance equation coefficient by coefficient through Kmax.
for k = 3:Kmax
    % Tune the central coefficients so the reduced leading term is -d*x^m.
    if tune && k <= m
        gS(k+1) = -sym(P.d)*(k == m) - a.'*HS(:,k-1);
    end
    acc = -B*HS(:,k-1) - w*(k == 3);
    for j = 2:k-2
        acc = acc + j*HS(:,j)*Fco(k+1-j, a, HS, gS);
    end
    HS(:,k) = Lam\acc;
end
for k = 2:m
    Fk = Fco(k, a, HS, gS);
    if k < m, assert(isAlways(Fk == 0), 'F_%d non nul', k); end
    if k == m, assert(isAlways(Fk == -sym(P.d)), 'F_m ~= -d'); end
end
end

function c = Fco(k, a, HS, gS)
c = a.'*HS(:,k-1);
if k+1 <= numel(gS), c = c + gS(k+1); end
end
