function f = general_proj_spec(IS_L_DATA,IS_FORWARD,G1_FLAG,in_matrix,sb,wa,num_rows,num_cols,num_bands,size_kernel)
% Last major changes:
% Normalize the kernel
% Linyun August 2016
% Try to remove the boundary issue
% Linyun 5/20/2016
tic;

    kernel_col = ones(num_cols, size_kernel);
    padding = ((size_kernel-1)/2);  
    cons_pre = in_matrix(:,:,1);
    cons_post = in_matrix(:,:,end);
    delta_tmp = mean(mean(diff(wa,1,2)));
    wa = [wa(:,1:padding) - delta_tmp*padding wa wa(:,end-padding+1:end) + delta_tmp*padding];
    
    % Values for band "area" obtained by (lambdaArray(final) - lambdaArray(initial)) / (num_bands - 1) for both L and S
    if IS_L_DATA
        band_area = 6.63363265;
    else
        band_area = 6.45642105;
    end
    
%     if IS_FORWARD
        fragment_col = ones(num_rows, num_cols, size_kernel);
        conv_sheet = ones(num_rows,num_cols,num_bands);
        num_iter = num_bands;
        fwhm_all = sb(1,:,:);
%     else
%         num_iter = num_bands + 2*padding;
%         factor = ones(num_rows,num_cols,num_bands);
%         conv_sheet_temp = zeros(num_rows, num_cols, num_iter+2*padding);
%         tmp = squeeze(sb(1,:,:));
%         delta_tmp = mean(mean(diff(tmp,1,2)));
%         fwhm_all = [tmp(:,1:padding) - delta_tmp*padding tmp tmp(:,end-padding+1:end) + delta_tmp*padding];
%     end
    
    for band = 1:num_iter
        if IS_L_DATA && ~G1_FLAG
            alpha1 = sb(2,:,band);
            gamma1 = sb(3,:,band);
            lambda1 = sb(4,:,band);
            alpha2 = sb(5,:,band);
            gamma2 = sb(6,:,band);
            lambda2 = sb(7,:,band);
            alpha3 = sb(8,:,band);
            gamma3 = sb(9,:,band);
            lambda3 = sb(10,:,band);
        else
            fwhm = fwhm_all(:,band)';
        end
        for band_iter = band-padding:band+padding
        %boundary condition
            if (band_iter<1)||(band_iter>num_iter)
%                 if IS_FORWARD
                    if band_iter < 1
                    fragment_col(:, :, band_iter-(band-padding)+1) = cons_pre;
                    else
                    fragment_col(:, :, band_iter-(band-padding)+1) = cons_post;
                    end
                    band_mark = band_iter + padding; 
%                 else
%                     kernel_col(:,band_iter-(band-padding)+1) = 0;

%                 end
            else
%                 if IS_FORWARD
                    band_mark = band_iter + padding;
                    %array of scene that kernel corresponds to
                    fragment_col(:, :, band_iter-(band-padding)+1) = in_matrix(:, :, band_iter);
%                 else
%                     %compute the same kernel as before
%                     %lambda_rep = repmat(lambda(band_iter), [1 num_cols]);
%                     lambda_rep = reshape(wa(:,band_iter),[1,num_cols]);
%                     if IS_L_DATA && ~G1_FLAG
%                         kernel_col(:,band_iter-(band-padding)+1) = band_area * (alpha1.*exp(-gamma1.*(lambda_rep-lambda1).^2) + alpha2.*exp(-gamma2.*(lambda_rep-lambda2).^2) + alpha3.*exp(-gamma3.*(lambda_rep-lambda3).^2) );
%                     else
%                         kernel_col(:, band_iter-(band-padding)+1) = band_area * double((2*sqrt(log(2))./(sqrt(pi).*fwhm)).*exp(-log(16)*((lambda_rep-wa(:,band)').^2)./(fwhm.^2)));
%                     end
%                     if band < padding + 1
%                         frag = cons_pre;
%                     elseif band >= padding + num_bands
%                         frag = cons_post;
%                     else
%                         frag = in_matrix(:,:,band-padding);
%                     end
%                 end             
            end
%             if IS_FORWARD
                if IS_L_DATA && ~G1_FLAG
                    kernel_col(:, band_iter-(band-padding)+1) = band_area * (alpha1.*exp(-gamma1.*(wa(:,band_mark)'-lambda1).^2) + alpha2.*exp(-gamma2.*(wa(:,band_mark)'-lambda2).^2) + alpha3.*exp(-gamma3.*(wa(:,band_mark)'-lambda3).^2) );
                else
                    kernel_col(:, band_iter-(band-padding)+1) = band_area * double((2*sqrt(log(2))./(sqrt(pi).*fwhm)).*exp(-log(16)*((wa(:,band_mark)'-wa(:,band+padding)').^2)./(fwhm.^2)));
                end
%             end
        end
%         if IS_FORWARD
            kernel_col = kernel_col./(repmat(sum(kernel_col,2),[1,size_kernel]));
            kernel = reshape(kernel_col, [1 num_cols size_kernel]);
            kernel_compact = repmat(kernel, [num_rows 1 1]);
            temp = dot_fast(fragment_col,kernel_compact,3);
            conv_sheet(:, :, band) = temp;
%         else
%             kernel = reshape(kernel_col, [1 num_cols size_kernel]);
%             kernel_compact = repmat(kernel, [num_rows 1 1]);
%             input = repmat(frag, [1 1 size_kernel]);
%             temp = input.*kernel_compact;
%             band_min_temp = band;
%             band_max_temp = band + 2*padding;
%             conv_sheet_temp(:, :,band_min_temp:band_max_temp) = conv_sheet_temp(:,:,band_min_temp:band_max_temp) + temp;
%         end
    end
    f = conv_sheet; 
    if IS_FORWARD
         %fill guess for measured scene with spectrally convolved data. Is this correct notation?
        display('Spectral Forward Projection: ');
    else
%         factor(:, :, :) = conv_sheet_temp(:,:, 1+2*padding:num_bands+2*padding);
%         'back_spec done.';   
%         f = factor; %output matrix
        display('Spectral Backprojection: ');
    end
    toc;
end