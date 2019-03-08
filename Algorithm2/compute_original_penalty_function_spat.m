
function [penalty_func] = compute_original_penalty_function_spat(sensitivity_image,new_estimate,beta_spat, delta_spat,  spat_penalty_on)        
% New version by Linyun 03/02/2016
% Speed up
    
    image_mask = sensitivity_image > 0;

    [c_num_rows,c_num_cols, num_bands] = size(new_estimate);
    penalty_function = zeros(c_num_rows,c_num_cols, num_bands);
    
        %first step, compute the spatial penalty
        if spat_penalty_on,
        
        new_estimate_padded = padarray(new_estimate, [1 1 0]);
        for row_inc =-1:1:1,
            for col_inc = -1:1:1,
            
                temp = new_estimate;
                neighbour_mask = image_mask;
                
                if row_inc == -1
                    temp = padarray(temp, [2 0 0], 'pre');
                    neighbour_mask = padarray(neighbour_mask, [2 0 0], 'pre');
                elseif row_inc == 1
                    temp = padarray(temp, [2 0 0], 'post');
                    neighbour_mask = padarray(neighbour_mask, [2 0 0], 'post');
                else
                    temp = padarray(temp, [1 0 0], 'both');
                    neighbour_mask = padarray(neighbour_mask, [1 0 0], 'both');
                end
            
                if col_inc == -1
                    temp = padarray(temp, [0 2 0], 'pre');
                    neighbour_mask = padarray(neighbour_mask, [0 2 0], 'pre');
                elseif col_inc == 1
                    temp = padarray(temp, [0 2 0], 'post');
                    neighbour_mask = padarray(neighbour_mask, [0 2 0], 'post');
                else
                    temp = padarray(temp, [0 1 0], 'both');
                    neighbour_mask = padarray(neighbour_mask, [0 1 0], 'both');
                end
            
                bracketterm = new_estimate_padded - temp;
                bracketterm = bracketterm(2:c_num_rows+1, 2:c_num_cols+1, :);%trim
                neighbour_mask = neighbour_mask(2:c_num_rows+1, 2:c_num_cols+1, :);
                neighbour_distance = sqrt(row_inc^2 + col_inc^2);
                if neighbour_distance>0,
                    neighbour_weight = 1/neighbour_distance;
                    
                    %this equation to compute the penalty function may not be correct        
                    temp2 = bracketterm/delta_spat;
                    temp = neighbour_weight*beta_spat*delta_spat^2.*log(cosh(temp2)).*neighbour_mask;
                    penalty_function = penalty_function + temp;
                end

            end
        end
        penalty_function(~image_mask) = 0;
        
        clearvars temp;
        clearvars new_estimate_padded;
        clearvars c_last_iter_padded;
        clearvars neighbour_mask;
        end
        
        
        
        penalty_func = penalty_function;
        
end