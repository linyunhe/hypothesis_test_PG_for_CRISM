function [penalty_func] = compute_decomposed_penalty_function_logcoshH(weight_spat,weight_spec,new_estimate,image_last_iter,beta_spat, delta_spat, beta_spec, delta_spec, spectral_span, spat_penalty_on, spec_penalty_on)        
% New version by Linyun 03/02/2016
% Speed up
    
    image_mask = weight_spat > 0;

    [c_num_rows,c_num_cols, num_bands] = size(image_last_iter);
    penalty_function = zeros(c_num_rows,c_num_cols, num_bands);
    
        %first step, compute the spatial penalty
        if spat_penalty_on,
        
        new_estimate_padded = padarray(new_estimate, [1 1 0]);
        c_last_iter_padded = padarray(image_last_iter, [1 1 0]);
        for row_inc =-1:1:1,
            for col_inc = -1:1:1,
            
                temp = image_last_iter;
                neighbour_mask = image_mask;
                temp_sen = weight_spat;
                
                if row_inc == -1
                    temp = padarray(temp, [2 0 0], 'pre');
                    temp_sen = padarray(temp_sen,[2 0 0],'pre');
                    neighbour_mask = padarray(neighbour_mask, [2 0 0], 'pre');
                elseif row_inc == 1
                    temp = padarray(temp, [2 0 0], 'post');
                    temp_sen = padarray(temp_sen,[2 0 0],'post');
                    neighbour_mask = padarray(neighbour_mask, [2 0 0], 'post');
                else
                    temp = padarray(temp, [1 0 0], 'both');
                    temp_sen = padarray(temp_sen, [1 0 0], 'both');
                    neighbour_mask = padarray(neighbour_mask, [1 0 0], 'both');
                end
            
                if col_inc == -1
                    temp = padarray(temp, [0 2 0], 'pre');
                    temp_sen = padarray(temp_sen, [0 2 0], 'pre');
                    neighbour_mask = padarray(neighbour_mask, [0 2 0], 'pre');
                elseif col_inc == 1
                    temp = padarray(temp, [0 2 0], 'post');
                    temp_sen = padarray(temp_sen, [0 2 0], 'post');
                    neighbour_mask = padarray(neighbour_mask, [0 2 0], 'post');
                else
                    temp = padarray(temp, [0 1 0], 'both');
                    temp_sen = padarray(temp_sen, [0 1 0], 'both');
                    neighbour_mask = padarray(neighbour_mask, [0 1 0], 'both');
                end
            
                bracketterm = 2*new_estimate_padded - c_last_iter_padded - temp;
                bracketterm = bracketterm(2:c_num_rows+1, 2:c_num_cols+1, :);%trim
                scalar = (weight_spat + temp_sen(2:c_num_rows+1,2:c_num_cols+1,:))/2;
                neighbour_mask = neighbour_mask(2:c_num_rows+1, 2:c_num_cols+1, :);
                neighbour_distance = sqrt(row_inc^2 + col_inc^2);
                if neighbour_distance>0,
                    neighbour_weight = 1/neighbour_distance;
                    
                    %this equation to compute the penalty function may not be correct        
                    temp2 = bracketterm/delta_spat;
                    temp = scalar.*neighbour_weight*beta_spat*delta_spat^2.*log(cosh(temp2)).*neighbour_mask;
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
        
        %second step, compute the spectral penalty
        if spec_penalty_on,
        for band_inc =-spectral_span:1:spectral_span,
        
            temp = image_last_iter;
            neighbour_mask = image_mask;
            temp_pnt = weight_spec;
            if band_inc<0,
                temp = padarray(temp, [0 0 -band_inc], 'replicate', 'pre');
                temp_pnt = padarray(temp_pnt, [0 0 -band_inc], 'replicate', 'pre');
                neighbour_mask = padarray(neighbour_mask, [0 0 -band_inc], 'pre');
            elseif band_inc>0,
                temp = padarray(temp, [0 0 band_inc], 'replicate', 'post');
                temp_pnt = padarray(temp_pnt, [0 0 band_inc], 'replicate', 'post');
                neighbour_mask = padarray(neighbour_mask, [0 0 band_inc], 'post');
            end

            bracketterm = 2*new_estimate - image_last_iter;
            scalar = weight_spec;
            
            if band_inc<0,
                bracketterm = bracketterm - temp(:,:,1:num_bands);
                scalar = scalar + temp_pnt(:,:,1:num_bands);
                neighbour_mask = neighbour_mask(:,:,1:num_bands);
            elseif band_inc>0,
                bracketterm = bracketterm - temp(:,:,1+band_inc:num_bands+band_inc);
                scalar = scalar + temp_pnt(:,:,1+band_inc:num_bands+band_inc);
                neighbour_mask = neighbour_mask(:,:,1+band_inc:num_bands+band_inc);
            end
            
            neighbour_distance = sqrt(band_inc^2);
            
            if neighbour_distance >0,
                neighbour_weight = 1/neighbour_distance;
            
                temp2 = bracketterm/delta_spec;
                temp = scalar/2.*neighbour_weight*beta_spec*delta_spec^2.*log(cosh(temp2)).*neighbour_mask;
                penalty_function = penalty_function + temp;
            
            end         
        end
        clearvars temp;
        clearvars temp2;
        clearvars neighbour_mask;
        end
        
        
        penalty_func = penalty_function;
        
end
