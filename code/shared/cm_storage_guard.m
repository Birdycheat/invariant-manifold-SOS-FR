function cm_storage_guard(pathname)
% Native MATLAB distribution: historical numerical algorithm, MIT project license.

% Resolve the path before enforcing the external-storage rule.
resolved = char(java.io.File(pathname).getCanonicalPath());
parts = regexp(lower(resolved), '[\\/]+', 'split');
assert(~any(startsWith(parts, 'onedrive')), 'centerjets:OneDriveOutput', ...
    'Use cm_solve_benchmark with an absolute output folder outside OneDrive.');
end
