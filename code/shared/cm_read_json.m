function value = cm_read_json(file)
%CM_READ_JSON Read JSON or gzipped JSON without changing reference files.
file = char(file);
if endsWith(file, '.gz', 'IgnoreCase', true)
    base = tempdir;
    assert(~contains(lower(base), 'onedrive'), 'cm:storage:OneDrive', ...
        'The system temporary directory must be outside OneDrive.');
    folder = tempname(base);
    mkdir(folder);
    guard = onCleanup(@() clearTemporary(folder));
    unpacked = gunzip(file, folder);
    assert(numel(unpacked) == 1 && startsWith(unpacked{1}, [folder, filesep]), ...
        'cm:json:archive', 'Unexpected gzip output.');
    value = jsondecode(fileread(unpacked{1}));
    clear guard;
else
    value = jsondecode(fileread(file));
end
end

function clearTemporary(folder)
assert(startsWith(folder, tempdir) && ~contains(lower(folder), 'onedrive'));
items = dir(folder);
for k = 1:numel(items)
    if ~items(k).isdir, delete(fullfile(folder, items(k).name)); end
end
rmdir(folder);
end
