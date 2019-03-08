function [o_penalty, d_penalty, new_image] = compute_logcosh_penalty_H_MSEE (sensitivity_image, spec_pnt,update_factor, image_last_iter, beta_spat, delta_spat, beta_spec, delta_spec, spectral_span, max_iter, spat_penalty_on, spec_penalty_on, use_mex)
tic;


    
    % 0 = H - update_factor*image_last_iter./c + weight*beta*delta*tanh{(2c-c(k)-c~(k))/delta}
    
    %check image_last_iter for negative values
%     if sum(sum(sum(image_last_iter<0)))
%         fprintf('Negative pixle values found in image_last_iter, please check input\n');
%         return;
%     end
%     
        % Run MEX version if we have it
    if use_mex
        if exist('compute_logcosh_penalty_H_mex','file')
            disp('Running C++ version of compute_logcosh_penalty_H_MSEE');
            global mex_num_threads;
            if(isempty(mex_num_threads))
                mex_num_threads=-1; % default to max threads
            end
            [o_penalty, d_penalty, new_image] = compute_logcosh_penalty_H_MSEE_mex(sensitivity_image, spec_pnt,update_factor, image_last_iter, beta_spat, delta_spat, beta_spec, delta_spec, spectral_span, max_iter, spat_penalty_on, spec_penalty_on, mex_num_threads);
            toc;
            return;
        else
            disp('* No binary for C++ compute_logcosh_penalty_H_mex present for this platform. Using Matlab version.');
        end
    else
        disp('Running Matlab version of compute_logcosh_penalty_H_MSEE');
    end
    [c_num_rows,c_num_cols, num_bands] = size(image_last_iter);
    image_mask = sensitivity_image > 0;
    
%     em_estimate_scaled = update_factor.*image_last_iter;
    %use the pure EM update as the initial value for the new estimate
    new_estimate = image_last_iter;
    
    if spat_penalty_on == 0 && spec_penalty_on == 0,
        %no penalty, pure EM update
        new_image = update_factor./sensitivity_image + image_last_iter;%no need to worry about NaN values in the result since they don't propogate during forward and backward projection
        d_penalty = 0;
        o_penalty = 0;
        return;
    end
    

    
    weight_spat = sensitivity_image;
    weight_spec = repmat(reshape(spec_pnt,[1,1,num_bands]),[c_num_rows,c_num_cols,1]).*sensitivity_image;
    
	h_val = 0.1;
	h_val_matrix = h_val*ones(c_num_rows,c_num_cols, num_bands);
    mask = update_factor < -10;

    iter = 1;

    while( iter<= max_iter),
        fprintf('Restricted Step Method, Iter: %d\n',iter);
        
        
        %force a nonnegativity constraint on the restricted region so it
        %ensures that the image update gets nonnegative value
        temp = new_estimate - h_val_matrix;
        h_val_matrix(temp<0) = 0.99*new_estimate(temp<0);
        
        
        
        penalty_function = zeros(c_num_rows,c_num_cols, num_bands);
        penalty_gradient = zeros(c_num_rows,c_num_cols, num_bands);
        penalty_hessian = zeros(c_num_rows,c_num_cols, num_bands);
        
        %first step, compute the spatial penalty
        if spat_penalty_on,
        new_estimate_padded = padarray(new_estimate, [1 1 0]);
        c_last_iter_padded = padarray(image_last_iter, [1 1 0]);
        for row_inc =-1:1:1,
            for col_inc = -1:1:1,
            
                temp = image_last_iter;
                temp_sen = weight_spat;
                neighbour_mask = image_mask;
            
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

                    temp2 = bracketterm/delta_spat;
                    
                    temp = scalar.*neighbour_weight*beta_spat*delta_spat^2.*log(cosh(temp2)).*neighbour_mask;
                    
                    penalty_function = penalty_function + temp;
                    
                    temp3 = exp(-2*temp2);
                    temp = scalar.*(2*neighbour_weight*beta_spat*delta_spat*(2./(1+temp3)-1).*neighbour_mask);
                    penalty_gradient = penalty_gradient + temp;
                    
                    temp = scalar.*(16*neighbour_weight*beta_spat*temp3./(1+temp3).^2.*neighbour_mask);
                    penalty_hessian = penalty_hessian + temp;
                end

            end
        end

        
        clearvars temp;
        clearvars temp2;
        clearvars temp3;
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
                temp = padarray(temp, [0 0 -band_inc],'replicate','pre');
                temp_pnt = padarray(temp_pnt, [0 0 -band_inc], 'replicate', 'pre');
                neighbour_mask = padarray(neighbour_mask, [0 0 -band_inc],  'pre');
            elseif band_inc>0,
                temp = padarray(temp, [0 0 band_inc], 'replicate', 'post');
                temp_pnt = padarray(temp_pnt, [0 0 band_inc],'replicate','post');
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
            
                temp3 = exp(-2*temp2);
                temp = scalar.*neighbour_weight.*beta_spec.*delta_spec.*(2./(1+temp3)-1).*neighbour_mask; 
                penalty_gradient = penalty_gradient + temp;
                
                temp = 8*scalar.*neighbour_weight.*beta_spec.*temp3./(1+temp3).^2.*neighbour_mask;
                penalty_hessian = penalty_hessian + temp;
            end         
        end
        clearvars temp;
        clearvars temp2;
        clearvars neighbour_mask;
        end

        
        penalty_function(~image_mask) = 0;
        penalty_gradient(~image_mask) = 0;
        penalty_hessian(~image_mask) = 0;
        
        tmp = image_last_iter - new_estimate;
        obj_function = sensitivity_image.*tmp.^2 + 2*update_factor.*tmp + penalty_function;
        obj_gradient = -2*sensitivity_image .* tmp - 2*update_factor + penalty_gradient;
        obj_heissian = 2*sensitivity_image + penalty_hessian;
        
        delta_raw = -obj_gradient./obj_heissian;    
        delta = delta_raw;
        hit_boundary = delta_raw>h_val_matrix | delta_raw<-h_val_matrix;
        delta(delta_raw>h_val_matrix) = h_val_matrix(delta_raw>h_val_matrix);
        delta(delta_raw<-h_val_matrix) = -h_val_matrix(delta_raw<-h_val_matrix);
        
        obj_function_quadratic_approximate = obj_function + obj_gradient.*delta + 0.5*obj_heissian.*delta.*delta;
        
        %evaluate the obj_function at new estimate and compute r
        rsm_new_estimate = new_estimate + delta;
        temp = image_mask.*compute_decomposed_penalty_function_logcoshH(weight_spat,weight_spec, rsm_new_estimate, image_last_iter,beta_spat, delta_spat, beta_spec, delta_spec,spectral_span, spat_penalty_on, spec_penalty_on);

        tmp = image_last_iter - rsm_new_estimate;
        obj_function_2 = sensitivity_image.*tmp.^2 + 2*update_factor.*tmp + temp;
        clearvars tmp;
        
        %compute r
        obj_diff = obj_function - obj_function_2;
        obj_approximate_diff = obj_function - obj_function_quadratic_approximate;
        r_vals = obj_diff./obj_approximate_diff;
        
        %if r_vals>0.75, then the quadratic approximates the function well.
        %and if the delta is within the restricted region, then the h_val is good enough
        
        %if the delta is at the boundary of the region, then we need to
        %increase the region
        h_val_matrix(r_vals>0.75 & hit_boundary) = 2*h_val_matrix(r_vals>0.75 & hit_boundary);
        
        
        %if r_vals<0.25, bad approximation, the curvature is too flat, we
        %need to use a smaller restricted region
        h_val_matrix(r_vals<0.25) = 0.25*abs(delta(r_vals<0.25));
        
        %else, h remains the same
        
        
        disp('Trust regions updated');
        
        
        mask = mask | r_vals > 0;
        %if r_vals>0, rsm_new_estimate accepted.
        new_estimate(r_vals>0) = rsm_new_estimate(r_vals>0);
        
        disp('Image estimate updated');
        
        iter = iter + 1;
        
        
        
    end%while(iter<max_iter)
    
    new_image = new_estimate;
    
    %compute the final penalty value
    temp = compute_decomposed_penalty_function_logcoshH(weight_spat,weight_spec, new_estimate, image_last_iter,beta_spat, delta_spat, beta_spec,delta_spec, spectral_span, spat_penalty_on, spec_penalty_on);
    d_penalty = sum(temp(~isnan(temp)));
    
    temp = weight_spat.*compute_original_penalty_function_spat(sensitivity_image, new_estimate, beta_spat, delta_spat,spat_penalty_on)+weight_spec.*compute_original_penalty_function_spec(sensitivity_image, new_estimate, beta_spec, delta_spec,spectral_span,spec_penalty_on);
    o_penalty = sum(temp(~isnan(temp)));
    
toc;
end
