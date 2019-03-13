function [mask,scalar,spatial_fig,scalar_fig] = cutpart(img,band,cutflag,wa_use)
   if length(band) > 3
       band = band(1:3);
   elseif length(band) < 3
      band = repmat(band(1),1,3); 
   end
   [num_row,num_col,num_band] = size(img);
   spatial_img = zeros(num_row,num_col,3);
%    spatial_img(:,:,1) = imadjust(img(:,:,band(1)));
%    spatial_img(:,:,2) = imadjust(img(:,:,band(2)));
%    spatial_img(:,:,3) = imadjust(img(:,:,band(3)));
   tmp = img(:,:,band(1));
   spatial_img(:,:,1) = imadjust(tmp/max(tmp(:)));
   tmp = img(:,:,band(2));
   spatial_img(:,:,2) = imadjust(tmp/max(tmp(:)));
   tmp = img(:,:,band(3));
   spatial_img(:,:,3) = imadjust(tmp/max(tmp(:)));
%    new = zeros(num_row,num_col);
   while cutflag == 1
       figure(1);
       imshow(spatial_img);
       title('Select a homogenous area by one click.');
       h = imfreehand(gca);
       new = createMask(h);
%    new = new | BW;
       mask = repmat(new,[1,1,num_band]);
       samples_num = sum(mask(:))/num_band;
       region = reshape(img(mask),[samples_num,num_band]);
       scalar_ra = zeros(1,num_band);
       for bands = 1:num_band
            tmp = region(:,bands);
            tmp = tmp(tmp > 0);
            mean_ra = mean(tmp);
            var_ra = var(tmp);
            scalar_ra(bands) = mean_ra/var_ra;
       end
    %
       scalar_fig = figure(2);
       plot(wa_use,1./scalar_ra);
       xlabel('Wavelength (\mu m)');
       ylabel('Scalar');
       title('Scalars for different bands');
       operate = inputdlg('Rechoose another area? Yes:Any number No:0');
       tmp = str2double(cell2mat(operate));
   
       
       if tmp ==0
           scalar = median(scalar_ra);
           spatial_fig = figure(1);
           spatial_img = show_withmask(img,mask,band,0.8);
           imshow(spatial_img);
           break;       
       end
   end
end