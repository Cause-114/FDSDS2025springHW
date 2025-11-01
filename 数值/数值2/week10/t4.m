% Define the function and its derivative and 2nd derivative.
f = @(x) exp(x) + log(x);
df = @(x) exp(x) + 1 ./ x;
d2f = @(x) exp(x) - 1 ./ x .^ 2;

% Set the range of x, the range of n, and the number of samples.
% And when n is equal to pln, plot the figure.
a = 1; b = 4; max_n = 20; pn = 200; pln = 4;
px = linspace(a, b, pn);
d1 = df(a); dn = df(b);

% Initialize the error of 1st and 2nd derivative.
err1 = zeros(max_n - 1, 1); err2 = zeros(max_n - 1, 1);
for n = 2:max_n
    x = linspace(a, b, n); y = f(x);
    [pdy, pd2y] = d1_spline(x, y, px, d1, dn);
    err1(n - 1) = get_max(df, px, pdy);
    err2(n - 1) = get_max(d2f, px, pd2y);
    if (n == pln)
        plot_d(df, px, pdy, 1);
        plot_d(d2f, px, pd2y, 2);
    end
end
plot_n(max_n, err1, err2);

% this fucntion is used to plot the original function 
% and the spline version of its derivative and their difference.
function plot_d(f, px, py, op)
    s = ["1st", "2nd"];
    rpy = f(px);
    figure;
    plot(px, rpy, px, py);
    title(s(op) + " Derivative of f(x)");
    xlabel("x"); ylabel(s(op) + " Derivative");
    legend("original", "spline version");
    figure;
    plot(px, abs(rpy - py));
    title("Error of " + s(op) + " Derivative");
    xlabel("x"); ylabel("Error");
end

% As its name, this function is used to get the maximum error 
% between the original function and the spline version.
function res = get_max(f, px, py)
    res = max(abs(f(px) - py));
end

% this fucntion is used to process and better pisplay 
% relationships between stepsize and max error of 1st and 2nd derivative.
function plot_n(n, err1, err2)
    figure; x = log10(3 ./ (1:n-1));
    err1 = log10(err1); err2 = log10(err2);
    plot(x, err1, x, err2);
    title("log10-version Error of 1st and 2nd Derivative");
    xlabel("log10-stepsize"); ylabel("log10-max Error");
    legend("1st Derivative", "2nd Derivative");
    % note here use the bound of 10\15, so your n must be larger than 15.
    fprintf("approximate slope of the first line:%.2f\n", (err1(10) - err1(15)) / (x(10) - x(15)));
    fprintf("approximate slope of the second line:%.2f\n", (err2(10) - err2(15)) / (x(10) - x(15)));
end


% this fucntion is used to get D1-cubic spline interpolation of f(x).
function [pdy, pd2y] = d1_spline(x, y, px, d1, dn)
    n = length(x);
    A = zeros(n - 1, 2);
    dy = zeros(n, 1);
    dy(1) = d1; dy(n) = dn; A(1, 1) = 1;
    s_l = (y(2) - y(1)) / (x(2) - x(1));
    % Initalize the A,b. And the upper triangular process.
    for i = 2:n - 1
        A(i, 2) = x(i) - x(i - 1); temp = x(i + 1) - x(i);
        A(i, 1) = 2 * (x(i + 1) - x(i - 1));
        s = (y(i + 1) - y(i)) / temp;
        dy(i) = 3 * (s_l * temp + s * A(i, 2));
        temp = -temp / A(i - 1, 1); A(i, 1) = A(i, 1) + A(i - 1, 2) * temp;
        dy(i) = dy(i) + temp * dy(i - 1); s_l = s;
    end
    % Solve the upper triangular system to get the derivatives.
    for i = n - 1:-1:2
        dy(i) = (dy(i) - A(i, 2) * dy(i + 1)) / A(i, 1);
    end

    id = 1; m = length(px);
    pdy = zeros(1, m); pd2y = zeros(1, m);
    for i = 1:n - 1
        h = x(i + 1) - x(i); s = (y(i + 1) - y(i)) / h;
        c = 3 * s / h - (2 * dy(i) + dy(i + 1)) / h;
        d = (dy(i + 1) + dy(i) - 2 * s) / (h ^ 2);
        while ((id < m && px(id) < x(i + 1)) || id == m)
            df = px(id) - x(i);
            pdy(id) = dy(i) + 2 * c * df + 3 * d * df ^ 2;
            pd2y(id) = 2 * c + 6 * d * df;
            id = id + 1;
        end
    end
end
