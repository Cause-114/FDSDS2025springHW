n = 20;
f1 = @(x) exp(x); f2 = @(x) x .^ 1.5;
res1 = exp(1) - 1; res2 = 0.4;
errf1 = zeros(1, n); errf2 = zeros(1, n);

for k = 1:n
    [x, w] = GLquadrature(k, 0, 1);
    approx1 = sum(w .* f1(x));
    approx2 = sum(w .* f2(x));
    errf1(k) = log10(abs(res1 - approx1));
    errf2(k) = log10(abs(res2 - approx2));
end

fprintf("log version of Gauss-Legendre");
fprintf("Quadrature Error in each turn of f1, f2\n");
disp(errf1)
disp(errf2)
plot(1:n, errf1, 1:n, errf2);
xlabel('Number of quadrature nodes');
ylabel('Absolute error (log scale)');
title('Gauss-Legendre Quadrature Error');
legend('f(x) = e^x', 'f(x) = x^{1.5}');
grid on;

function [x, w] = GLquadrature(N, a, b)
    beta = 1 ./ sqrt(4 - (1:N - 1) .^ (-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x = diag(D);
    [x, i] = sort(x);
    x = 0.5 * (b - a) * x + 0.5 * (b + a);
    w = ((b - a) * V(1, i) .^ 2)';
end
