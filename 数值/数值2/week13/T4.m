% Solve the PDE:
% u_t = u_xx, 0 < x < 1, t > 0
% u(0, t) = 1 - 2t, u(1, t) = -2t
% u(x, 0) = 1 - x^2

EndT = 1;
% Here, it's better to make ht can divide EndT
hx = 0.1; ht = 0.001;
nx = round(1 / hx) + 1;
nt = round(EndT / ht) + 1;
X = linspace(0, 1, nx);
T = linspace(0, EndT, nt);
% 4 different methods to solve this "initial value problem".
Explict_Euler(X, T);
Implict_Euler(X, T);
Trapzoidal(X, T);
RK4(X, T);

function Explict_Euler(X, T)
    nx = length(X); hx = X(2) - X(1);
    nt = length(T); ht = T(2) - T(1);
    Y = zeros(nx, nt);
    Y(:, 1) = 1 - X .^ 2;
    for i = 2:nt
        Y(1, i) = 1 - 2 * T(i); Y(nx, i) = -2 * T(i);
        d2x = get_2derivative(Y(2:nx - 1, i - 1), Y(1, i - 1), Y(nx, i - 1), hx);
        Y(2:nx - 1, i) = Y(2:nx - 1, i - 1) + ht * d2x;
    end
    figure; imagesc(T, X, Y);
    xlabel('t'); ylabel('x');
    colorbar; title(sprintf("Explict Euler Method, with hx = %.2f, ht = %.3f", hx, ht));
end

function Implict_Euler(X, T)
    nx = length(X); hx = X(2) - X(1);
    nt = length(T); ht = T(2) - T(1);
    Y = zeros(nx, nt);
    Y(:, 1) = 1 - X .^ 2;
    for i = 2:nt
        Y(1, i) = 1 - 2 * T(i); Y(nx, i) = -2 * T(i);
        A = spdiags([-1, 2 + hx ^ 2 / ht, -1], [-1, 0, 1], nx - 2, nx - 2);
        b = hx ^ 2 / ht * Y(2:nx - 1, i - 1); b(1) = b(1) + Y(1, i); b(end) = b(end) + Y(nx, i);
        Y(2:nx - 1, i) = A \ b;
    end
    figure; imagesc(T, X, Y);
    xlabel('t'); ylabel('x');
    colorbar; title(sprintf("Implict Euler Method, with hx = %.2f, ht = %.3f", hx, ht));
end

function Trapzoidal(X, T)
    nx = length(X); hx = X(2) - X(1);
    nt = length(T); ht = T(2) - T(1);
    Y = zeros(nx, nt);
    Y(:, 1) = 1 - X .^ 2;
    for i = 2:nt
        Y(1, i) = 1 - 2 * T(i); Y(nx, i) = -2 * T(i);
        A = spdiags([-1, 2 + 2 * hx ^ 2 / ht, -1], [-1, 0, 1], nx - 2, nx - 2);
        b = 2 * hx ^ 2 / ht * Y(2:nx - 1, i - 1); b(1) = b(1) + Y(1, i); b(end) = b(end) + Y(nx, i);
        for j = 2:nx - 1
            b(j - 1) = b(j - 1) + Y(j + 1, i - 1) - 2 * Y(j, i - 1) + Y(j - 1, i - 1);
        end
        Y(2:nx - 1, i) = A \ b;
    end
    figure; imagesc(T, X, Y);
    xlabel('t'); ylabel('x');
    colorbar; title(sprintf("Trapzoidal Method, with hx = %.2f, ht = %.3f", hx, ht));
end

function RK4(X, T)
    nx = length(X); hx = X(2) - X(1);
    nt = length(T); ht = T(2) - T(1);
    Y = zeros(nx, nt);
    Y(:, 1) = 1 - X .^ 2;
    for i = 2:nt
        Y(1, i) = 1 - 2 * T(i); Y(nx, i) = -2 * T(i);
        tmp = Y(2:nx - 1, i - 1); y0 = Y(1, i - 1); yn = Y(nx, i - 1);
        k1 = ht * get_2derivative(tmp, y0, yn, hx);
        k2 = ht * get_2derivative(tmp + k1 / 2, y0 - ht, yn - ht, hx);
        k3 = ht * get_2derivative(tmp + k2 / 2, y0 - ht, yn - ht, hx);
        k4 = ht * get_2derivative(tmp + k3, y0 - 2 * ht, yn -2 * ht, hx);
        Y(2:nx - 1, i) = tmp + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
    figure; imagesc(T, X, Y);
    xlabel('t'); ylabel('x');
    colorbar; title(sprintf("RK4 Method, with hx = %.2f, ht = %.3f", hx, ht));
end

function y = get_2derivative(y, y0, yn, hx)
    y = [y0; y(:); yn];
    nx = length(y);
    for i = 2:nx - 1
        y(i - 1) = (y(i + 1) - 2 * y(i) + y(i - 1)) / (hx ^ 2);
    end
    y = y(1:nx - 2);
end
