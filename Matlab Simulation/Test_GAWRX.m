n = 10;
N = 2^n;
SNR = 3;
sigma = sqrt(1/10^(SNR/10));
x_mean = 2/(sigma^2) * ones(N,1);

u_mean = zeros(N,1); 
u_mean = GA_Construction_punctured_polar_codes(x_mean,N,n);
[Z_s, Z_index] = sort(u_mean);

index = [1:N]';
% [Z_s, Z_index, index]
Z_s(N/2+1)