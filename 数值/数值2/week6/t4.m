% The problem mentioned "Try to fit the data using a periodic cubic 
% spline with equispaced nodes {1,4,7,...,19,22,25}.".It's easy to use 
% interpolation to find for periodic cubic spline. But how to find the
% so-called "periodic cubic spline" by fitting? In my perspective, 
% this problem needs us to use the least square method to find the 
% coefficients of the periodic cubic spline. And note that the periodic 
% cubic spline is completely determined by the value at each node, so we 
% can savely just use the f-value as the parameter for Guass-Newton iteration.
x = [1, 3, 5, 7, 8, 9, 10, 11, 12, 13, ...
         14, 15, 16, 17, 18, 19, 20, 21, 22, 23];
y = [36.37, 36.23, 36.21, 36.26, 36.38, 36.49, 36.60, 36.63, 36.66, 36.68, ...
         36.69, 36.73, 36.74, 36.78, 36.82, 36.84, 36.87, 36.86, 36.77, 36.59];

px = linspace(0, 48, 500);
res = period_fit(x, y, px); figure;
plot(x, y, 'o', px, res, '-');
legend('data', 'fit result')
xlabel('Time / h', 'Interpreter', 'latex');
ylabel("temperature / $\,^{\circ}\mathrm{C}$", "Interpreter", "latex");
title('Periodic Spline fit for temperature data');

x(end + 1) = x(1) + 24; y(end + 1) = y(1);
res = period_interpolate(x, y, px); figure;
plot(x(1:end - 1), y(1:end - 1), 'o', px, res, '-');
legend('data', 'interploate result')
xlabel('Time / h', 'Interpreter', 'latex');
ylabel("temperature / $\,^{\circ}\mathrm{C}$", "Interpreter", "latex");
title('Periodic spline interpolation for temperature data');

function res = period_fit(x, y, px)
    iy = [36.37; 36.21; 36.26; 36.6; 36.68; 36.74; 36.84; 36.69];
    n = length(iy); m = length(x);
    A = 4 * eye(n, n); A(1, n) = 1; A(n, 1) = 1;
    for i = 1:n - 1
        A(i, i + 1) = 1; A(i + 1, i) = 1;
    end
    A = inv(A); tmp0 = A(:, 1); tmp1 = A(:, n);
    for i = 1:n - 1
        tmp = A(:, mod(i - 2, n) + 1) - A(:, i + 1);
        A(:, mod(i - 2, n) + 1) = tmp1; tmp1 = tmp;
    end
    A(:, n) = A(:, n - 1) - tmp0; A(:, n - 1) = tmp1;
    % the deraivatives of each node can be obtained by A * iy.
    J = zeros(m, n);
    for i = 1:m
        id = floor((x(i) + 2) / 3); d = x(i) - 3 * id + 2;
        c1 = 1 -5/27 * d ^ 2 +2/27 * d ^ 3; c2 = d -4/9 * d ^ 2 +1/9 * d ^ 3;
        c3 = 5/27 * d ^ 2 -2/27 * d ^ 3; c4 = (d ^ 3 - d ^ 2) / 9;
        J(i, :) = c2 * A(id, :) + c4 * A(mod(id, n) + 1, :);
        J(i, id) = J(i, id) + c1; J(i, mod(id, n) + 1) = J(i, mod(id, n) + 1) + c3;
    end
    iy = J \ y'; dy = A * iy;
    iy(end + 1) = iy(1); dy(end + 1) = dy(1);
    [px, ord] = sort(mod(px - 1, 24) + 1);
    id = 1; m = length(px); res = zeros(1, m);
    for i = 1:n
        s = (iy(i + 1) - iy(i)) / 3;
        c = 3 * s / 3 - (2 * dy(i) + dy(i + 1)) / 3;
        d = (dy(i + 1) + dy(i) - 2 * s) / 9;
        while ((id < m && px(id) < i * 3 + 1) || id == m)
            diff = px(id) - i * 3 + 2;
            res(ord(id)) = iy(i) + dy(i) * diff + c * diff ^ 2 + d * diff ^ 3;
            id = id + 1;
        end
    end
end

function res = period_interpolate(x, y, px)
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
    for i = 1:n - 1
        h = x(i + 1) - x(i); s = (y(i + 1) - y(i)) / h;
        c = 3 * s / h - (2 * dy(i) + dy(i + 1)) / h;
        d = (dy(i + 1) + dy(i) - 2 * s) / (h ^ 2);
        while ((id < m && px(id) < x(i + 1)) || id == m)
            diff = px(id) - x(i);
            res(ord(id)) = y(i) + dy(i) * diff + c * diff ^ 2 + d * diff ^ 3;
            id = id + 1;
        end
    end
end