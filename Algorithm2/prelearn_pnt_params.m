function [f, spectrogram] = prelearn_pnt_params(data,N_row,N_col,length_window,scalar)
   tmp = multiple_pixel_spectral(data,N_row,N_col,length_window);
   spectrogram = tmp;
   tmp = tmp(end,:);
   tmp = (tmp-min(tmp))/(max(tmp)-min(tmp));
   tmp = (normcdf(tmp,0,median(tmp))-0.5)*(scalar-1)*2+1;  
   d1 = designfilt('lowpassiir','FilterOrder',12, ...
       'HalfPowerFrequency',0.15,'DesignMethod','butter');
   f = filtfilt(d1,tmp);
end

