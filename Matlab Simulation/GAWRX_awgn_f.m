 function [free_sets,frozen_sets,punctured_sets] = GAWRX(n,Nr,K,SNR)
%% 参数设置
% 设定仿真参数
N = 2^n;
N_remain = Nr;
R = K/N_remain;
N_punct = N - N_remain;
sigma = sqrt(1/10^(SNR/10));

N_temp = zeros(N,1);
N_temp = rvsl([1:N]');
puncturing_positions = N_temp(N - N_punct + 1:N); %模拟打孔的位置
puncturing_positions = sort(puncturing_positions);
%% 计算P
x_mean = zeros(N,1);
x_mean(:,1) = 2/(sigma^2);
x_mean(puncturing_positions,1) = 300;%被打孔位置可靠性无穷大
%信息位的均值
u_mean = zeros(N,1); 
u_mean = GA_Construction_punctured_polar_codes(x_mean,N,n);
u_sigma2 = 2*u_mean;
pe_u =1 - qfunc(-u_mean./(sqrt(u_sigma2)+0.0000000001));
pe_u(N-N_punct+1:N) = 0.5;
Z_temp_p = pe_u;
[Z_s, Z_index] = sort(Z_temp_p);
free_sets=Z_index(1:K);
free_sets=sort(free_sets);
frozen_sets=Z_index(K+1:N);
frozen_sets=sort(frozen_sets);
punctured_sets=puncturing_positions;
end
 