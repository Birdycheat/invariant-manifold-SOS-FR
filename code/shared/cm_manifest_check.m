function result = cm_manifest_check(root)
%CM_MANIFEST_CHECK Verify every declared file and reject undeclared files.
file = fullfile(root,'MANIFEST.json');
assert(isfile(file),'cm:Manifest','MANIFEST.json is missing.');
manifest = jsondecode(fileread(file));
expected = string({manifest.files.bundle_path});
% Check declared bytes and reject unlisted files as well as missing ones.
for k=1:numel(manifest.files)
    row = manifest.files(k);
    target = fullfile(root,strrep(row.bundle_path,'/',filesep));
    assert(isfile(target) && strcmp(cm_sha256(target),row.sha256), ...
        'cm:Integrity','Missing or changed reference file: %s.',row.bundle_path);
end
actual = walk(root,'');
assert(isempty(setxor(sort(actual),sort([expected,"MANIFEST.json"]))), ...
    'cm:Integrity','Undeclared files are present in the reference repository.');
result = struct('files',numel(expected),'manifest_sha256',cm_sha256(file));
end

function names = walk(root, relative)
names = strings(1,0);
items = dir(fullfile(root,relative));
for k=1:numel(items)
    name = items(k).name;
    if any(strcmp(name,{'.','..'})) || (isempty(relative) && strcmp(name,'.git')), continue; end
    part = fullfile(relative,name);
    if items(k).isdir
        names = [names,walk(root,part)]; %#ok<AGROW>
    else
        names(end+1) = string(strrep(part,filesep,'/')); %#ok<AGROW>
    end
end
end
