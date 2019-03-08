function [ isOkay, message ] = has_mex( )
%HAS_MEX Determines if CRISM has valid, usable binaries on this system
%   Returns true only if MEX binaries are present for this system (judged
%   by trying to run mex_check, which is made in CRISM MLM's MexMakeFile) &
%   they can actually run here (i.e., no missing dependencies or anything).

if ~exist('mex_check','file')
    isOkay = false;
    message = 'Missing compiled binaries compatible with current operating system';
    
    return;
end

try
    % run a MEX C++ function compilied just like ones we really need
    res = mex_check();
catch e
    isOkay = false;
    message = sprintf('Trying to run MEX function got error %s',e.message);
    return;
end


isOkay = true;
message = '';

end

