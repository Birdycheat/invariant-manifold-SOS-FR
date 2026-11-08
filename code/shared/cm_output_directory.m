function folder = cm_output_directory(root, requested)
%CM_OUTPUT_DIRECTORY Keep generated files outside the reference and OneDrive.
if strlength(string(requested))==0, requested = tempname; end
requested = char(requested);
assert(java.io.File(requested).isAbsolute(), 'cm:OutputPath', 'OutputDir must be absolute.');
% Resolve aliases before checking the repository and OneDrive boundaries.
folder = char(java.io.File(requested).getCanonicalPath());
root = char(java.io.File(root).getCanonicalPath());
parts = regexp(lower(folder), '[\\/]+', 'split');
assert(~any(startsWith(parts,'onedrive')), 'cm:OneDriveOutput', 'OutputDir must be outside OneDrive.');
assert(~strcmpi(folder,root) && ~startsWith(lower(folder),[lower(root) filesep]), ...
    'cm:ReferenceOutput', 'OutputDir must be outside the reference repository.');
if ~isfolder(folder), [ok,msg] = mkdir(folder); assert(ok,'cm:OutputPath','%s',msg); end
end
