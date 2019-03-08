function [ B ] = safeSubset( A, row_min, row_max, col_min, col_max, defVal)
%SAFESUBSET Subsets a 3D matrix with default values for invalid indices
%   Copies a portion of a 3D matrix A into the matrix B. If all indices for A are valid,
%   this function behaves the same as B = A(row_min:row_max,col_min:col_max,:);
%   Cells with out-of-bounds rows and columns are instead replaced by the scalar given
%   as defVal.

B = defVal * ones(row_max-row_min+1,col_max-col_min+1,size(A,3));

min_viable_row = max([1 row_min]);
max_viable_row = min([row_max size(A,1)]);

min_viable_col = max([1 col_min]);
max_viable_col = min([col_max size(A,2)]);

B_min_insert_row = min_viable_row-row_min+1;
B_max_insert_row = B_min_insert_row + (max_viable_row-min_viable_row);

B_min_insert_col = min_viable_col-col_min+1;
B_max_insert_col = B_min_insert_col + (max_viable_col-min_viable_col);

B(B_min_insert_row:B_max_insert_row,B_min_insert_col:B_max_insert_col,:) = A(min_viable_row:max_viable_row,min_viable_col:max_viable_col,:);

end
