function  f = general_proj_spat_precomputed(huge_kernel, IS_FORWARD, in_matrix, num_rows,num_cols,num_bands,size_kernel,lat_number,lon_number,c_num_rows,c_num_cols)
% Computes the kernel. Forward projection if IS_FORWARD is true, back projection otherwise.
% New version by Linyun 03/02/2016 passing huge_kernel as a variable.

    tic;
    padding = ((size_kernel-1)/2); 
   
    if IS_FORWARD
        conv_layer = ones(num_rows,num_cols,num_bands);
        fragment = ones(size_kernel,size_kernel,num_bands);
        num_elements = size_kernel*size_kernel;
    else
        conv_layer = zeros(c_num_rows+2*padding,c_num_cols+2*padding,num_bands);
        factor = ones(c_num_rows,c_num_cols,num_bands);    
    end
    
    %Compute the kernel for the spatial TF
    for col = 1:num_cols
        for row = 1:num_rows
            if IS_FORWARD
                % build fragment from elements within row_iter and col_iter ranges. These were previously defined as:
                % row_iter = lat_number(row,col)-padding:lat_number(row,col)+padding
                % col_iter = lon_number(row,col)-padding:lon_number(row,col)+padding  
                % This method fills the rest of the space with zeros.
                fragment = safeSubset(in_matrix,lat_number(row,col)-padding,lat_number(row,col)+padding,lon_number(row,col)-padding,lon_number(row,col)+padding,0);
            end
            
            total_kernel = huge_kernel(:,:,:,row,col);
            if IS_FORWARD
                kernel_temp = reshape(total_kernel, num_elements, num_bands);
                fragment_temp = reshape(fragment, num_elements, num_bands);
                conv_layer(row,col,:) = sum(kernel_temp.*fragment_temp,1); % dot product
            else
                temp2 = repmat (in_matrix(row, col, 1:num_bands), [size_kernel size_kernel 1]);
                temp = temp2.*total_kernel;
                %boundaries of patch in the new image that the kernel
                %corresponds to
                row_min = lat_number(row,col);
                row_max = lat_number(row,col)+2*padding;
                col_min = lon_number(row,col);
                col_max = lon_number(row,col)+2*padding;

                %compute section of scene
                conv_layer(row_min:row_max,col_min:col_max, :) = conv_layer(row_min:row_max,col_min:col_max, :) + temp;       
            end 
        end
    end
    if IS_FORWARD
        f = conv_layer;    
        display('Spatial Forward Projection: ');
    else
        factor(:, :, :) = conv_layer(1+padding:c_num_rows+padding, 1+padding:c_num_cols+padding, :);
        'back_spat done.';
        f = factor;
        display('Spatial Backprojection: ');
    end
    toc;
end

