N = 10;
n = log2(N);
bin = zeros(1,n);
pe_u = zeros(1,N);
for i = 1:N
    ii=i;
    for j=1:n
        index = 2^(n-j);
        bin(j) = floor((ii-1)/index);
        ii = ii-bin(j)*index;
    end
    
    for j=1:n
        pe_u(i)=pe_u(i)+bin(n+1-j)*2^((j-1)/4);
    end
end
Z_temp_p = pe_u;
[Z_s, Z_index] = sort(Z_temp_p);
Z_index'