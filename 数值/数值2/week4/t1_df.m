% This code is used to calculate the high order derivative of f 
% by using the divided differences method, in order to estimate the 
% upper bound of the error.
n = 20;
f = @(x)(1 + 25 * x ^ 2) ^ -1;
x = linspace(-0.1, 0.1, n + 1);
y = arrayfun(f, x);

for i = 1:n
    tmp2 = y(i);
    for j = 1:n - i + 1
        tmp1 = (y(j + i) - y(i + j - 1)) / (x(j + i) - x(j));
        y(j + i - 1) = tmp2;
        tmp2 = tmp1;
    end
    y(n + 1) = tmp2;
end

figure
plot(0:n, log10(abs(y)))
title('High Order Derivative of f near 0(approximated by differentiation)')
xlabel('Order of Derivative')
ylabel('$log_{10}(|f^{(order)}|)$', 'Interpreter', 'latex')