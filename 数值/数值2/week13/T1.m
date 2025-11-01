FigureMaker1(10);
FigureMaker2(round(logspace(1,3,100)));

function r = f_actual(x)
    % Actual function value
    e = exp(1); a = 2 * e / (1 + e); b = 2 - a;
    r = x .^ 2 + 2 - a * exp(-x) - b * exp(x);
end

function b = finite_difference(n)
    % Solve the ODE: -u''(x)+u(x)=x^2 by finite difference method on the
    % interval [0,1] with u(0)=0,u(1)=1; n: number of small intervals;
    h = 1 / n;
    b = (0:1 / n:1) .^ 2;
    a1 = 2 / h ^ 2 + 1; a2 = -1 / h ^ 2;
    b(n) = b(n) - a2;
    A = ones(n - 1, 2) * a2; A(1, 1) = a1; A(n - 1, 2) = 0;
    % Because A is a tridiagonal matrix, we can solve it by upper-triangulize it.
    for i = 2:n - 1
        tmp = -a2 / A(i - 1, 1);
        A(i, 1) = a1 + A(i - 1, 2) * tmp;
        b(i + 1) = b(i + 1) + b(i) * tmp;
    end
    for i = n - 1:-1:1
        b(i + 1) = (b(i + 1) - A(i, 2) * b(i + 2)) / A(i, 1);
    end
end

function b = finite_element(n)
    % Solve the ODE: -u''(x)+u(x)=x^2 by finite element method on the
    % interval [0,1] with u(0)=0,u(1)=1; n: number of small intervals;
    h = 1 / n;
    b = ((0:n) .^ 2 +1/6) .* h ^ 3;
    b(1) = 0; b(n + 1) = 1;
    a1 = 2 * (n + h / 3); a2 = -n + h / 6;
    b(n) = b(n) - a2;
    A = ones(n - 1, 2) * a2; A(1, 1) = a1; A(n - 1, 2) = 0;
    % Because A is a tridiagonal matrix, we can solve it by upper-triangulize it.
    for i = 2:n - 1
        tmp = -a2 / A(i - 1, 1);
        A(i, 1) = a1 + A(i - 1, 2) * tmp;
        b(i + 1) = b(i + 1) + b(i) * tmp;
    end
    for i = n - 1:-1:1
        b(i + 1) = (b(i + 1) - A(i, 2) * b(i + 2)) / A(i, 1);
    end
end

function FigureMaker1(n)
    x = linspace(0, 1, n + 1);
    r = f_actual(x);
    r1 = finite_difference(n);
    r2 = finite_element(n);
    % function value figure
    figure;
    plot(x, r1, x, r2, x, r);
    title(sprintf("Function Value, where n=%d", n));
    xlabel('x'); ylabel('u(x)'); grid on;
    legend('Finite Difference', 'Finite Element', 'Actual');
    % difference figure
    figure;
    plot(x, abs(r1 - r), 'r--', x, abs(r2 - r), 'g--');
    title(sprintf("Difference, where n=%d", n));
    xlabel('x'); ylabel('\Delta u(x)'); grid on;
    legend('FD Error', 'FE Error');
end

function FigureMaker2(n_list)
    n = length(n_list);
    residual1 = zeros(n, 1);
    residual2 = zeros(n, 1);
    for i = 1:n
        n_ = n_list(i);
        x = linspace(0, 1, n_ + 1);
        r = f_actual(x);
        r1 = finite_difference(n_);
        r2 = finite_element(n_);
        residual1(i) = log(norm(r1 - r, inf));
        residual2(i) = log(norm(r2 - r, inf));
    end
    % log-log scale of relationship between difference and n
    n_list = log(n_list);
    figure;
    plot(n_list, residual1, n_list, residual2);
    title("log-log scale of relationship between difference and n");
    xlabel('log(n)'); ylabel('log(max(\Delta u))'); grid on;
    legend('Finite Difference', 'Finite Element');
    % least square fit of the two line
    [a1, b1] = least_squares(n_list, residual1);
    [a2, b2] = least_squares(n_list, residual2);
    fprintf("Below is the least square fit of the two lines in FigureMaker2\n");
    fprintf('FD line : y = %.2f x + %.2f\n', a1, b1);
    fprintf('FM line : y = %.2f x + %.2f\n', a2, b2);
end

function [a, b] = least_squares(x, y)
    n = length(x);
    x_ = mean(x); y_ = mean(y);
    a = (dot(x, y) - n * x_ * y_) / (dot(x, x) - n * x_ * x_);
    b = (y_ - a * x_);
end
