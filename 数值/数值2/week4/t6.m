x1 = [-0.93, -0.74, -0.67, -0.50, -0.34, -0.26, ...
          -0.18, -0.11, -0.08, -0.06, -0.03, 0.02, 0.13, 0.16, ...
          0.20, 0.21, 0.26, 0.30, 0.40, 0.43, 0.44, 0.47, 0.48, ...
          0.53, 0.69, 0.75];
y1 = [0.13, 0.24, -0.12, 0.26, 1.27, 0.46, ...
          0.14, 0.00, 0.15, 0.56, 0.94, 1.45, 1.31, 0.36, ...
          0.00, 0.23, 0.71, 1.23, 1.45, 0.64, 0.33, 0.03, -0.09, ...
          0.29, 0.7, 0.31];
x2 = [-0.93, -0.80, -0.72, -0.51, -0.31, 0.00, 0.32, 0.48, 0.54, 0.66, 0.75];
y2 = [0.13, -0.70, -0.95, -1.24, -1.44, -1.56, -1.23, -0.88, -0.69, -0.36, 0.31];
px = linspace(-0.93, 0.75, 1000);

res = d2_spline(x1, y1, px);
plot(px, res, '--')
hold on
res = d2_spline(x2, y2, px);
plot(px, res)
hold on
plot(x1, y1, 'o', x2, y2, 'o');
title("hand-shape-spline")

function res = d2_spline(x, y, px)
    n = length(x);
    A = zeros(n - 1, 2);
    dy = zeros(n, 1);
    A(1, 1) = 1; A(1, 2) = 2;
    s_l = (y(2) - y(1)) / (x(2) - x(1));
    dy(1) = 3 * s_l;
    % Initalize the A,b. And the upper triangular process.
    for i = 2:n - 1
        A(i, 2) = x(i) - x(i - 1); temp = x(i + 1) - x(i);
        A(i, 1) = 2 * (x(i + 1) - x(i - 1));
        s = (y(i + 1) - y(i)) / temp;
        dy(i) = 3 * (s_l * temp + s * A(i, 2));
        temp = -temp / A(i - 1, 1); A(i, 1) = A(i, 1) + A(i - 1, 2) * temp;
        dy(i) = dy(i) + temp * dy(i - 1); s_l = s;
    end
    dy(n) = (3 * s_l - (dy(n - 1)) / A(n - 1, 1)) / (2 - A(n - 1, 2) / A(n - 1, 1));
    % Solve the upper triangular system to get the dy.
    for i = n - 1:-1:1
        dy(i) = (dy(i) - A(i, 2) * dy(i + 1)) / A(i, 1);
    end
    id = 1; m = length(px); res = zeros(1, m);
    for i = 1:n - 1
        h = x(i + 1) - x(i); s = (y(i + 1) - y(i)) / h;
        c = 3 * s / h - (2 * dy(i) + dy(i + 1)) / h;
        d = (dy(i + 1) + dy(i) - 2 * s) / (h ^ 2);
        while ((id < m && px(id) < x(i + 1)) || id == m)
            diff = px(id) - x(i);
            res(id) = y(i) + dy(i) * diff + c * diff ^ 2 + d * diff ^ 3;
            id = id + 1;
        end
    end
end
