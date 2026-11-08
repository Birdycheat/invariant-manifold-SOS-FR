function write_center(fn, instance, coords, enc, rstr, c, composite)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.
% Preserve matrix dimensions and column-major order for rational reconstruction.
flat = @(A) struct('dims', size(A), 'vals', {arrayfun(@(q) sprintf('%.17g', q), A(:)', 'uni', 0)});
J = struct('instance', instance, 'coords', coords, 'enc', enc, 'r', rstr, 'composite', composite, ...
           'ok', c.ok, 't', c.t, 'e', sprintf('%.17g', c.e), 'iter', c.iter, 'time', c.time, ...
           'E1', c.E{1}, 'E2', c.E{2}, 'E3', c.E{3}, 'Q1', flat(c.Q1));
J.Q2 = cellfun(flat, c.Q2(:)', 'uni', 0);
J.Q3 = cellfun(flat, c.Q3(:)', 'uni', 0);
fid = fopen(fn, 'w'); fprintf(fid, '%s', jsonencode(J)); fclose(fid);
end
