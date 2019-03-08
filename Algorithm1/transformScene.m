function [ transformedScene] = transformScene( c, varargin )
%TRANSFORMSCENE Rotate Scene for saving to file
%   Rotates a CRISM scene 180 degrees so that rows go from North to South
%   and columns go from West to East. This is done to make display in ENVI
%   work as expected (i.e., North upwards)
%
%   OPTIONS:
%   - If 'inverse' is given as an argument after the scene, the inverse of
%   the transformation is applied instead. While at this time the
%   transformation (180 deg rotation) is its own inverse, this futureproofs
%   against needing to change the code in multiple places later if this
%   changes.

% If any argument in varargin is the string 'inverse', run the inverse
% transformation instead. Ignore other args.

is_inverse = ~isempty(varargin) && any(strcmp('inverse',varargin));

if is_inverse % tranform from saving -> working orientation
    
    transformedScene = flip(flip(c,1),2); % flip over both spatial axes
    
else % transform from working -> saving orientation
    
    transformedScene = flip(flip(c,1),2); % flip over both spatial axes
    
end

end

