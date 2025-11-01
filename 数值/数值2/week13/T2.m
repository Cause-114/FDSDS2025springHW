% We can solve the generalized egvalue problem by computing B^-1*A(here B is 
% triangular, so O(n^2) time is taken). At the same time, maybe we can use 
% some "sophisticated" method to just obtain the m-th smallest eigvalue.
% But I'm a 'lazy' guy, so I just ask "eigs" for help.

% Left parameter stands for the number of samples right one stands for number 
% of eigvalues you want to obtain. It's the same for all the functions in this file;
FigureMaker1(100, 3);

function [v, f] = finite_difference(n, m)
    % Consider -u''(x)=\lambda u(x) on [0, pi] with u(0)=0,u(pi)=0; computing some small 
    % \lambda makes the solution not equal to zero all the time by finite difference method.
    h = pi / n;
    A = spdiags([-1 / h ^ 2, 2 / h ^ 2, -1 / h ^ 2], [-1, 0, 1], n - 1, n - 1);
    [f, v] = eigs(A, m, 'smallestabs');
    [v, s] = sort(diag(v));
    f = f(:, s);
    % regularize the f to make the biggest value equal to one, the same bellow. 
    for i = 1:m
        id = ceil(n / 2 / m);
        if f(id, i) >= 0
            f(:, i) = f(:, i) / norm(f(:, i), inf);
        else
            f(:, i) = -f(:, i) / norm(f(:, i), inf);
        end
    end
    f = [zeros(1, m); f; zeros(1, m)];
end

function [v, f] = finite_element(n, m)
    % Consider -u''(x)=\lambda u(x) on [0, pi] with u(0)=0,u(pi)=0; computing some small
    % \lambda makes the solution not equal to zero all the time by finite element method.
    h = pi / n;
    A = spdiags([-1 / h, 2 / h, -1 / h], [-1, 0, 1], n - 1, n - 1);
    B = spdiags([h / 6, 2 * h / 3, h / 6], [-1, 0, 1], n - 1, n - 1);
    [f, v] = eigs(A, B, m, 'smallestabs');
    [v, s] = sort(diag(v));
    f = f(:, s);
    for i = 1:m
        id = ceil(n / 2 / m);
        if f(id, i) >= 0
            f(:, i) = f(:, i) / norm(f(:, i), inf);
        else
            f(:, i) = -f(:, i) / norm(f(:, i), inf);
        end
    end
    f = [zeros(1, m); f; zeros(1, m)];
end

function FigureMaker1(n, m)
    x = linspace(0, pi, n + 1);
    [v1, f1] = finite_difference(n, m);
    [v2, f2] = finite_element(n, m);
    for i = 1:m
        r = sin(i * x);
        % function value figure
        figure;
        plot(x, f1(:, i)', x, f2(:, i)', x, r);
        title(sprintf("Function Value(n=%d ev1: %.1f, ev2: %.1f, actual: %d)", n, v1(i), v2(i), i ^ 2));
        xlabel('x'); ylabel('u(x)'); grid on;
        legend('Finite Difference', 'Finite Element', 'Actual');
        % difference figure
        figure;
        plot(x, abs(f1(:, i)' - r), '--', x, abs(f2(:, i)' - r), '.');
        title(sprintf("Difference(n=%d ev1: %.1f, ev2: %.1f, actual: %d)", n, v1(i), v2(i), i ^ 2));
        xlabel('x'); ylabel('\Delta u(x)'); grid on;
        legend('FD Error', 'FE Error');
    end
end
