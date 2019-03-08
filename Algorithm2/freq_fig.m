function f = freq_fig(data,L)

N = length(data);
% L = 71;
span_half = (L-1)/2;
w = hamming(L);
num_part = N;
f = zeros(span_half+2,num_part);  
for i = 1:num_part
    if (i-span_half)<1
        d = [data(span_half-i+2:-1:2);data(1:1:i+span_half)];
    elseif (i+span_half)>N
        d = [data(i-span_half:1:N);data(N-1:-1:2*N-i-span_half)];
    else
        d = [data(i-span_half:1:i+span_half)];
    end
    d = d.*w; 
    y = abs(fft(d,L+2)).^2;
    f(:,i) = y(1:1:span_half +2 )/max(y);
end

end