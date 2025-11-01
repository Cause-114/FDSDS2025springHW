n = 100;
f = @(x, y) exp(x) * sin(y);
err_list = zeros(1, n); res = 0;
tans = (exp(1) + cos(1) - sin(1)) / 2 - 1;
for i = 1:100
    res = main(0, 1, 0, 1, f, i);
    err_list(i) = log10(abs(res - tans));
end
plot((1:n).^0.5, err_list);
xlabel('Number of quadrature nodes');
ylabel('Absolute error (log scale)');
title('composite quadrature rule over 2d triangle');
grid on;
fprintf("Integral of f(x, y) = %.10g\n", res);

function res = main(x1, x2, y1, y2, f, n)
    res = 0; hx = (x2 - x1) / n; hy = (y2 - y1) / n;
    for i = 1:n
        x0 = x1; y0 = y1 + (i - 1) * hy;
        for j = 1:i - 1
            res = res + cat(x0, x0 + hx, y0, y0 + hy, f);
            res = res + cat(x0 + hx, x0, y0, y0 - hy, f);
            x0 = x0 + hx; y0 = y0 - hy;
        end
        res = res + cat(x0, x0 + hx, y0, y0 + hy, f);
    end
end

function y = cat(x1, x2, y1, y2, f)
    S = 1/2 * (x2 - x1) * (y2 - y1);
    tmp = 5/6 * x1 +1/6 * x2;
    x2 = 1/3 * x1 +2/3 * x2; x1 = tmp;
    tmp = 5/6 * y1 +1/6 * y2;
    y2 = 1/3 * y1 +2/3 * y2; y1 = tmp;
    y = (f(x1, y1) + f(x2, y1) + f(x1, y2)) / 3 * S;
end
