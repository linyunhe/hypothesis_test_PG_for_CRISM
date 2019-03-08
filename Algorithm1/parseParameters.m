function [ argStruct ] = parseParameters( arglist )
%PARSEPARAMETERS Parse a list of parameters into a structure
%   Given a cell array of even length composed of parameter, value pairs,
%   returns a structure with fields for each parameter and its value. This
%   is easier to parse with common Matlab functions than looping through
%   the parameters every time.

if mod( size(arglist), 2 ) ~= 0
    error('requires even number of arguments (in order param, value, ..., param, value, etc');
end

argStruct = cell2struct(arglist(2:2:end),arglist(1:2:end),2);

end

