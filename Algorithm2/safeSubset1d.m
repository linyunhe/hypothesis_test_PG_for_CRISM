function [ B ] = safeSubset1d( A, pos_min, pos_max, defVal)
%SAFESUBSET1D Subsets a 1D matrix with default values for invalid indices
%   Copies a portion of a 1D matrix A into the matrix B. If all indices for A are valid,
%   this function behaves the same as B = A(row_min:row_max);.
%   Cells with out-of-bounds rows and columns are instead replaced by the scalar given
%   as defVal.

% keep at size 1 whichever dimension of the input has size 1 (support both
% row and column vectors)
if size(A,1)==1
    B = defVal * ones(1,pos_max-pos_min+1);
elseif size(A,2)==1
    B = defVal * ones(pos_max-pos_min+1,1);
else
    error('Input matrix is not 1D');
end

min_viable_pos = max([1 pos_min]);
max_viable_pos = min([pos_max length(A)]);

B_min_insert_pos = min_viable_pos-pos_min+1;
B_max_insert_pos = B_min_insert_pos + (max_viable_pos-min_viable_pos);

B(B_min_insert_pos:B_max_insert_pos) = A(min_viable_pos:max_viable_pos);

end

