function [d_kl,scalar] = Calculate_best_scalar(G,scalar_min,scalar_max)
   
   while abs(scalar_min-scalar_max)/scalar_max > 0.05
       display(sprintf('Searching from %f to %f',scalar_min,scalar_max));
       scalar_mean = (scalar_min + scalar_max)/2;
       G1 = G * scalar_min;
       pvalue = 1 - chi2cdf(G1(:),1);
       [N,~] = hist(pvalue,100);
       pdf = N/sum(N);
       d_kl_min = sum(pdf.*log(pdf)) + log(100);
       G1 = G * scalar_max;
       pvalue = 1 - chi2cdf(G1(:),1);
       [N,~] = hist(pvalue,100);
       pdf = N/sum(N);
       d_kl_max = sum(pdf.*log(pdf)) + log(100);
       G1 = G * scalar_mean;
       pvalue = 1 - chi2cdf(G1(:),1);
       [N,~] = hist(pvalue,100);
       pdf = N/sum(N);
       d_kl_mean = sum(pdf.*log(pdf)) + log(100);
       d_kl_real_min = min([d_kl_min,d_kl_max,d_kl_mean]);
       if d_kl_min == d_kl_real_min
            scalar_max = scalar_min;
            scalar_min = scalar_max - (scalar_mean - scalar_max);
       elseif d_kl_max == d_kl_real_min
            scalar_min = scalar_max;
            scalar_max = scalar_min + (scalar_min - scalar_mean);
       else
            scalar_max = (scalar_mean+scalar_max)/2;
            scalar_min = (scalar_mean+scalar_min)/2;
       end
        
   end
   scalar = (scalar_min + scalar_max)/2;
   G1 = G * scalar;
   pvalue = 1 - chi2cdf(G1(:),1);
   [N,~] = hist(pvalue,100);
   pdf = N/sum(N);
   d_kl = sum(pdf.*log(pdf)) + log(100);
   display('Searching Done!')
end