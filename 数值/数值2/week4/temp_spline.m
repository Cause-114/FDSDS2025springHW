% I think it's better to make the period 24h instead of 23h, so I add one more
% point at 25h, thus the f, df, d^2f at 25h are the same
% as at 1h. At the same time, this change make best use of the data points, by
% avoiding lose the data points at 23h or 1h to make sure they are the same.
x = [1, 3, 5, 7, 8, 9, 10, 11, 12, 13, ...
         14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25];
y = [36.37, 36.23, 36.21, 36.26, 36.38, 36.49, 36.60, 36.63, 36.66, 36.68, ...
         36.69, 36.73, 36.74, 36.78, 36.82, 36.84, 36.87, 36.86, 36.77, 36.59, 36.37];

px = linspace(0, 48, 500);
res = period_spline(x, y, px);
plot(x, y, 'o', px, res, '-');
legend('data', 'spline result')
xlabel('Time / h', 'Interpreter', 'latex');
ylabel("temperature / $\,^{\circ}\mathrm{C}$", "Interpreter", "latex");
title('Periodic Spline for  temperature data');

function res = period_spline(x, y, px)
    % x, y: sample points, y(i) = f(x(i))
    % px: the point for figure plot
    % res: the result of spline interpolation at px
    n = length(x);
    A = eye(n - 1, n - 1);
    b = zeros(n - 1, 1);
    s_l = (y(n) - y(n - 1)) / (x(n) - x(n - 1));
    h_l = x(n) - x(n - 1);
    for i = 1:n - 1
        s = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
        l = mod(i - 2, n - 1) + 1; r = mod(i, n - 1) + 1;
        A(i, r) = h_l; A(i, l) = x(i + 1) - x(i);
        A(i, i) = 2 * (h_l + A(i, l));
        b(i) = 3 * (s_l * A(i, l) + s * h_l);
        s_l = s; h_l = A(i, l);
    end
    dy = A \ b; dy(end + 1) = dy(1);
    [px, ord] = sort(mod(px - x(1), x(n) - x(1)) + x(1));
    id = 1; m = length(px); res = zeros(1, m);
    h = x(n) - x(n - 1); s = (y(n) - y(n - 1)) / h;
    c = 3 * s / h - (2 * dy(n - 1) + dy(n)) / h;
    d = (dy(n - 1) + dy(n) - 2 * s) / (h ^ 2);
    d1f = dy(n - 1) + 2 * c * h + 3 * d * h ^ 2;
    d2f = 6 * d * h + 2 * c;
    for i = 1:n - 1
        h = x(i + 1) - x(i); s = (y(i + 1) - y(i)) / h;
        c = 3 * s / h - (2 * dy(i) + dy(i + 1)) / h;
        d = (dy(i + 1) + dy(i) - 2 * s) / (h ^ 2);
        while ((id < m && px(id) < x(i + 1)) || id == m)
            diff = px(id) - x(i);
            res(ord(id)) = y(i) + dy(i) * diff + c * diff ^ 2 + d * diff ^ 3;
            id = id + 1;
        end
        % Bellow is the condition check part, I'm quite sure my code above
        % is right for that I've tested it under many cases.
        % If you find it looks ugly, feel free to delete it.
        if (abs(dy(i) - d1f) > 1e-10)
            fprintf('Warning: The spline is not C^1 continuous at i=%d\n', i);
            fprintf("say last df=%f, while this one is %f\n", d1f, dy(i))
        end
        if (abs(2 * c - d2f) > 1e-10)
            fprintf('Warning: The spline is not C^2 continuous at i=%d\n', i);
            fprintf("say last d2f=%f, while this one is %f\n", d2f, 2 * c)
        end
        d1f = dy(i) + 2 * c * h + 3 * d * h ^ 2;
        d2f = 6 * d * h + 2 * c;
    end
end
