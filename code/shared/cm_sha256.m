function digest = cm_sha256(file)
%CM_SHA256 Hash a file without changing its bytes.
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(file,'rb');
assert(fid>=0,'cm:Input','Cannot read %s.',file);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    bytes = fread(fid,1024*1024,'*uint8');
    if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
end
bytes = typecast(md.digest(),'uint8');
digest = lower(reshape(dec2hex(bytes,2).',1,[]));
end
