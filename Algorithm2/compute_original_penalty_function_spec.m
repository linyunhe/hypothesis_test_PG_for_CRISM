
function [penalty_func] = compute_original_penalty_function_spec(sensitivity_image,new_estimate,beta_spec, delta_spec, spectral_span, spec_penalty_on)        
% New version by Linyun 03/02/2016
% Speed up
    
    image_mask = sensitivity_image > 0;

    [c_num_rows,c_num_cols, num_bands] = size(new_estimate);
    penalty_function = zeros(c_num_rows,c_num_cols, num_bands);
    
        if spec_penalty_on,
        for band_inc =-spectral_span:1:spectral_span,
        
            temp = new_estimate;
            neighbour_mask = image_mask;
            
            if band_inc<0,
                temp = padarray(temp, [0 0 -band_inc], 'replicate','pre');
                neighbour_mask = padarray(neighbour_mask, [0 0 -band_inc], 'pre');
            elseif band_inc>0,
                temp = padarray(temp, [0 0 band_inc],  'replicate','post');
                neighbour_mask = padarray(neighbour_mask, [0 0 band_inc],  'post');
            end

            bracketterm = new_estimate;
            
            if band_inc<0,
                bracketterm = bracketterm - temp(:,:,1:num_bands);
                neighbour_mask = neighbour_mask(:,:,1:num_bands);
            elseif band_inc>0,
                bracketterm = bracketterm - temp(:,:,1+band_inc:num_bands+band_inc);
                neighbour_mask = neighbour_mask(:,:,1+band_inc:num_bands+band_inc);
            end
            
            neighbour_distance = sqrt(band_inc^2);
            
            if neighbour_distance >0,
                neighbour_weight = 1/neighbour_distance;
            
                %this equation to compute the penalty function may not be correct            
                temp2 = bracketterm/delta_spec;
                temp = neighbour_weight*beta_spec*delta_spec^2.*log(cosh(temp2)).*neighbour_mask; 
                 
                penalty_function = penalty_function + temp;            
            end         
        end
        penalty_function(~image_mask) = 0;
        
        clearvars temp;
        clearvars temp2;
        clearvars neighbour_mask;
        end
        
        
        
        penalty_func = penalty_function;
        
end