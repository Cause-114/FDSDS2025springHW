f = @(x) (x - 2) ^ 9;
g = @(x) x ^ 9 - 18 * x ^ 8 + 144 * x ^ 7 - 672 * x ^ 6 + ...
    2016 * x ^ 5 - 4032 * x ^ 4 + 5376 * x ^ 3 - 4608 * x ^ 2 + 2304 * x - 512;
x = linspace(1.9, 2.1, 1000);
y = zeros(length(x), 1);
z = zeros(length(x), 1);

for i = 1:length(x)
    y(i) = f(x(i));
    z(i) = g(x(i));
end
plot(x, y);
title("function displayed not in the expanded form")
figure;
plot(x, z);
title("function displayed in the expanded form")
