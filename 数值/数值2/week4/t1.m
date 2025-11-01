% Actually, I'm not sure what's the exact meaning of 'hermite spline'
% in the problem statement, but I guess it means that the interplote
% function has the same f and 1st derivative values at the interpolation points,
% so called "piecewise interpolation", not the "cubic spline",
% which makes sure that the 1st and 2nd derivative values of interpolated function are continuous.

n = 7;% feel free to change the number of interpolation points here.
f = @(x)(1 + 25 * x ^ 2) ^ -1;
df = @(x)(-50 * x) * (1 + 25 * x ^ 2) ^ (-2);
plot_x = linspace(-1, 1, 2001);
plot(plot_x, arrayfun(f, plot_x)); hold on;

x = linspace(-1, 1, n);
y = arrayfun(f, x);
res = interplot_Newton(y, x, plot_x);
plot(plot_x, res); hold on;

dy = arrayfun(df, x);
res = interplot_Hermit(y, dy, x, plot_x);
plot(plot_x, res); hold on;

x_cheb = cos(pi / (2 * n):pi / n:pi);
x_cheb = sort(x_cheb);
y_cheb = arrayfun(f, x_cheb);
res_cheb = interplot_Newton(y_cheb, x_cheb, plot_x);
plot(plot_x, res_cheb, '--'); hold on;
legend('Original', 'equispaced polynomial', 'Hermit spline', 'Chebyshev polynomial');
title(sprintf("number of interpolation points = %d", n))


function res = interplot_Newton(f, x, plot_x)
    % f(i) stores the function value at the i-th interpolation point.
    % x stand for the list of interpolation points.(i.e. x_1, x_2, ..., x_n)
    % plot_x is the x-axis of the plot, and res is the interpolated function value at plot_x.
    n = length(x);
    for i = 1:n - 1
        % each turn, we caculate the i+1-th order differential
        % and at this time, only the front n-i part stores valid data.
        tmp2 = f(i);
        for j = 1:n - i
            tmp1 = (f(j + i) - f(i + j - 1)) / (x(j + i) - x(j));
            f(j + i - 1) = tmp2;
            tmp2 = tmp1;
        end
        f(n) = tmp2;
    end
    res = zeros(length(plot_x), 1);
    for i = 1:length(plot_x)
        prod_temp = 1;
        for j = 1:n - 1
            res(i) = res(i) + f(j) * prod_temp;
            prod_temp = prod_temp * (plot_x(i) - x(j));
        end
        res(i) = res(i) + f(n) * prod_temp;
    end
end

function res = interplot_Hermit(f, df, x, plot_x)
    % f is the function value at x, df is the derivative value at x.
    % x is the interpolation point, plot_x is the x-axis of the plot,
    % and res is the interpolated function value at plot_x.
    n = length(x); m = length(plot_x);
    res = zeros(m, 1);
    id = 1; c = zeros(4, 1);
    for i = 1:n - 1
        % bellow caculate the coefficients of the cubic hermite interpolation polynomial.
        % similar to the caulation of differential in Newton interpolation.
        c(1) = f(i); c(2) = df(i); c(4) = df(i + 1);
        c(3) = (f(i + 1) - f(i)) / (x(i + 1) - x(i));
        c(4) = (c(4) - c(3)) / (x(i + 1) - x(i));
        c(3) = (c(3) - c(2)) / (x(i + 1) - x(i));
        c(4) = (c(4) - c(3)) / (x(i + 1) - x(i));
        % c(1) = f(i); c(2) = df(i);
        % c(3) = ((f(i + 1) - f(i)) / (x(i + 1) - x(i)) - df(i)) / (x(i + 1) - x(i));
        % c(4) = ((df(i + 1) + df(i)) * (x(i + 1) - x(i)) - 2 * (f(i + 1) - f(i))) / (x(i + 1) - x(i)) ^ 3;
        while ((id < m && plot_x(id) < x(i + 1)) || (id == m))
            res(id) = c(1) + c(2) * (plot_x(id) - x(i));
            res(id) = res(id) + c(3) * (plot_x(id) - x(i)) ^ 2;
            res(id) = res(id) + c(4) * (plot_x(id) - x(i)) ^ 2 * (plot_x(id) - x(i + 1));
            id = id + 1;
        end
    end
end
