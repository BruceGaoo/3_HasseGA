x = 2:1:15;
y = zeros(size(x));
z = zeros(size(x));
for i=1:size(x,2)
    y(i) = floor((6/pi)^0.5*2^x(i)/x(i)^1.5);
    z(i) = ceil(log2(y(i)));
end