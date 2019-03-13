function f = show_withmask(img,mask,band,density)
% show overlap with red color
  if size(img)~=size(mask)
     error('Inputs must have the same size');
  end
  if length(band) ~= 3
     error('Please choose three bands as RGB'); 
  end
  
  [num_row, num_col, num_band] = size(img);
  showimg = zeros(num_row,num_col,3);
  tmp = img(:,:,band(1));
  showimg(:,:,1) = imadjust(tmp/max(tmp(:)) + density * mask(:,:,band(1)));
  tmp = img(:,:,band(2));
  showimg(:,:,2) = imadjust(tmp/max(tmp(:)));
  tmp = img(:,:,band(3));
  showimg(:,:,3) = imadjust(tmp/max(tmp(:)));
  f = showimg;
end