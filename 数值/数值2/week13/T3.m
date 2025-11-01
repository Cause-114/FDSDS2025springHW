% According to the large scale of the problem(100*100 order of linear 
% equations), maybe it's a better choice to use inexact method like CG.
% But I'm a 'lazy' guy, so I just ask "\" for help.
n = 100;
U = finite_differ_2d(n);
xi = linspace(-1, 1, n + 1);
[X, Y] = meshgrid(xi, xi);
figure; surf(X, Y, U); xlabel('x'); ylabel('y');
title("3D-vision of the function");
test(u);

function test(U)
    n = size(U, 1) - 1;
    L = zeros(n - 1, n - 1);
    for i = 1:n - 1
        for j = 1:n - 1
            temp = U(i, j + 1) + U(i + 2, j + 1) + U(i + 1, j) + U(i + 1, j + 2);
            L(i, j) = (4 * U(i + 1, j + 1) - temp) * n ^ 2;
        end
    end
    % figure; imagesc(L); colorbar; title("discrete Laplace operater on each points");
    if (norm(L, inf) > 1e-6 * n)
        fprintf("Value of Laplace operator is not zero at each points\n");
        disp(norm(L, inf))
    else
        fprintf("under 1e-6 tolerance, no problem!\n");
    end
end

function U = finite_differ_2d(n)
    F = zeros(n - 1, n - 1);
    F(1, :) = (2 / n:2 / n:2 - 2 / n);
    F(n - 1, :) = F(1, :);
    F(:, 1) = F(:, 1) + ((1 - 2 / n:-2 / n:2 / n - 1) .^ 2 - 1)';
    F(:, n - 1) = F(:, n - 1) + ((1 - 2 / n:-2 / n:2 / n - 1) .^ 2 + 1)';
    T = 2 * eye(n - 1);
    for i = 1:n - 2
        T(i + 1, i) = -1;
        T(i, i + 1) = -1;
    end
    T = kron(eye(n - 1), T) + kron(T', eye(n - 1));
    F = F(:);
    u = T \ F;
    U = zeros(n + 1, n + 1);
    U(1, 1:n + 1) = (0:2 / n:2);
    U(n + 1, 1:n + 1) = U(1, 1:n + 1);
    U(2:n, 1) = ((1 - 2 / n:-2 / n:2 / n - 1) .^ 2 - 1)';
    U(2:n, n + 1) = U(2:n, 1) + 2;
    U(2:n, 2:n) = reshape(u, n - 1, n - 1);
end
