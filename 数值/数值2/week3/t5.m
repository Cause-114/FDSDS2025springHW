xi = [1.00000, 0.80902, 0.30902, -0.30902, -0.80902, -1.00000, -0.80902, -0.30902, 0.30902, 0.80902];
yi = [0.00000, 0.58779, 0.95106, 0.95106, 0.58779, 0.00000, -0.58779, -0.95106, -0.95106, -0.58779];
zi = [-1.00000, -2.6807, 5.6161, 5.6161, -2.6807, -1.00000, -2.6807, 5.6161, 5.6161, -2.6807];
point = complex(xi, yi);
res = interplot_Newton(zi, point);
x = linspace(-1, 1, 100);
y = linspace(-1, 1, 100);
z = zeros(100, 100);
for i = 1:100
    for j = 1:100
        z(i, j) = res(complex(x(j), y(i)));
    end
end

[X, Y] = meshgrid(x, y);
% fillter out the points not in the unit disk.
mask = (X .^ 2 + Y .^ 2) <= 1;
z(~mask) = NaN;
figure;
imagesc(x, y, real(z));
colorbar;
title('the interpolated function using Newton method');
xlabel('x');
ylabel('y');
figure;
surf(X, Y, real(z));
title('3d surface plot of the interpolated function using Newton method');
xlabel('x');
ylabel('y');

% the interpolated function used in T4.
function res = interplot_Newton(f_list, x_list)
    % f_list(i) stores the function value at the i-th interpolation point. 
    % x_list stand for the list of interpolation points.(i.e. x_1, x_2, ..., x_n)
    n = length(x_list);
    a_list = zeros(n, 1); a_list(1) = f_list(1);
    % a_list(i) stores f[x_1, ..., x_i]
    for i = 1:n - 1
    % each turn, we caculate the i+1-th order differential
    % and at this time, only the front n-i part stores valid data.
        temp = (f_list(2) - f_list(1)) / (x_list(1 + i) - x_list(1));
        a_list(i + 1) = temp;
        for j = 2:n - i
            f_list(j - 1) = temp;
            temp = (f_list(j + 1) - f_list(j)) / (x_list(j + i) - x_list(j));
        end
        f_list(n - i) = temp;
    end
    % the same as: $\sum_{i=0}^{n-1} (f[x_1, ..., x_{i+1}]\prod_{j=1}^{i} (t - x_j))$
    res = @(t) sum(arrayfun(@(i) a_list(i) * prod(t - x_list(1:i - 1)), 1:n));
end
