% 28 Haz 2006 Erdal Arikan
% Encoder for the folding construction
% This function assumes that u is of dimension N=2^n, in other words u has
% data bits as well as frozen bits which may have been set arbitrarily.

function x = encode(u);
N = size(u,1); % N must be a power of 2
n = log2(N);
if n==1
    x = [mod(u(1)+u(2),2); u(2)];
    return;
else
    x1 = encode(mod(u(1:N/2,1)+u(N/2+1:N,1),2));
    x2 = encode(u(N/2+1:N,1));
    x = [x1; x2];
end
%x=rvsl(x); % bit-reversal This must be applied externally
return
    