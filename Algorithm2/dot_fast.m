function [ result ] = dot_fast( a, b, dim )
%DOT_FAST Takes the dot-product of two vectors/matrices without frills.
%   Like the builtin dot(), but without the all the safety checks. Does not
%   attempt to support complex numbers or give nice-looking error messages
%   like dot() does (instead, the multiplication operation itself will
%   yell).

if nargin==2,
  result = sum(a.*b);
else
  result = sum(a.*b,dim);
end

end

