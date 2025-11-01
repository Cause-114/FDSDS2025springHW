px = linspace(0, pi * 2, 1000);
plot(px, sin(px), '--', 'DisplayName', 'sin(x)')
hold on;

n = 6; % number of sample points, of course, you can change it.
x = linspace(0, pi * 2, n); y = sin(x);
dn = 1; % change the boundary derivatives conditions d1 and dn as you like.
for d1 = linspace(-1, 1, 5)
    res = d1_spline(x, y, d1, dn, px);
    name = sprintf(" spline with f^'(0)=%.2f", d1);
    plot(px, res, 'DisplayName', name);
    hold on;
end
plot(x, sin(x), 'o', 'DisplayName', 'sample points');
title("spline results of sin(x) with different $f^{'}(0)$", 'Interpreter', 'latex')
legend show;

function res = d1_spline(x, y, d1, dn, px)
    % x, y: sample points, y(i) = f(x(i))
    % d1, dn: boundary derivatives conditions
    % px: the point for figure plot
    % res: the result of spline interpolation at px
    % make sure the x and px is in ascending order, otherwise 
    % you should add the sort part.
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
    % Solve the upper triangular system to get the dy.
    for i = n - 1:-1:2
        dy(i) = (dy(i) - A(i, 2) * dy(i + 1)) / A(i, 1);
    end

    id = 1; m = length(px); res = zeros(1, m);
    for i = 1:n - 1
        h = x(i + 1) - x(i);s = (y(i + 1) - y(i)) / h; 
        c = 3 * s / h - (2 * dy(i) + dy(i + 1)) / h;
        d = (dy(i + 1) + dy(i) - 2 * s) / (h ^ 2);
        while ((id < m && px(id) < x(i + 1)) || id == m)
            diff = px(id) - x(i);
            res(id) = y(i) + dy(i) * diff + c * diff ^ 2 + d * diff ^ 3;
            id = id + 1;
        end
        % Bellow is the condition check part, I'm quite sure my code above
        % is right for that I've tested it under many cases. For that I haven't 
        % changed the bound conditions d1 and dn, so I don't check it.
        % If you find it looks ugly, feel free to delete it. 
        if (i > 1 && abs(dy(i) - d1f) > 1e-10)
            fprintf('Warning: The spline is not C^1 continuous at i=%d\n', i);
            fprintf("say last df=%f, while this one is %f\n", d1f, dy(i))
        end
        if (i > 1 && abs(2 * c - d2f) > 1e-10)
            fprintf('Warning: The spline is not C^2 continuous at i=%d\n', i);
            fprintf("say last d2f=%f, while this one is %f\n", d2f, 2 * c)
        end
        d1f = dy(i) + 2 * c * h + 3 * d * h ^ 2;
        d2f = 6 * d * h + 2 * c;
    end
end

% A = eye(n, n);
% b = zeros(n, 1);
% b(1) = d1; b(n) = dn;
% s_l = (y(2) - y(1)) / (x(2) - x(1));

% for i = 2:n - 1
%     s = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
%     A(i, i + 1) = x(i) - x(i - 1); A(i, i - 1) = x(i + 1) - x(i);
%     A(i, i) = 2 * (x(i + 1) - x(i - 1));
%     b(i) = 3 * (s_l * A(i, i - 1) + s * h_l);
%     s_l = s;
% end

% dy = A \ b;
