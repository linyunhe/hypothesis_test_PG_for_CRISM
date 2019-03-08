function f = multiple_pixel_spectral(data,N_row,N_col,length_window)
display('Pecomputing penalty parameters...');
tic_pre = tic;

mask = data > 0;
mask = mean(mask,3)>0;
[num_row,num_col,num_bands] = size(data);
spec_eng_red = zeros((length_window-1)/2+2,num_bands,10);

N = 100;
spec_eng_mean = zeros((length_window-1)/2+2,num_bands,10);
for m = 1:N
    k=1;
    rand_row = ceil(rand(1,N_row)*num_row);
    rand_col = ceil(rand(1,N_col)*num_col);
    for i = 1:N_row
        for j = 1:N_col
            if mask(rand_row(i),rand_col(j))
                spec_eng_red(:,:,k) =  freq_fig(squeeze(data(rand_row(i),rand_col(j),:)),length_window);
                k = k+1;
            end
        end
        spec_eng_mean(:,:,m) = mean(spec_eng_red,3);
    end
end
f = mean(spec_eng_mean,3);

display('Pecomputing penalty parameters...done');
toc(tic_pre);
end