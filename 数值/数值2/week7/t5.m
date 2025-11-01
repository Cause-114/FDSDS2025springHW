f = @(x) exp(x);
x = [-0.8, -0.3, 0.4, 0.7]; tol = 1e-6;
[x1, c1, err1] = remez(f, f, x, -1, 1, tol);
c1(1) = c1(1) - c1(2) * x1(1) + c1(3) * x1(1) * x1(2);
c1(2) = c1(2) - c1(3) * (x1(1) + x1(2));
g1 = @(x) c1(3) * x .^ 2 + c1(2) * x + c1(1);

[x2, c2, err2] = cross_ver(tol);
g2 = @(x) c2(3) * x .^ 2 + c2(2) * x + c2(1);

fprintf("Remez algorithm:\n")
fprintf("crossing point set: %.6f, %.6f, %.6f, %.6f\n", x1(1), x1(2), x1(3), x1(4));
fprintf("fitting polynomial: %.6fx^2+%.6fx+%.6f\n", c1(3), c1(2), c1(1));
fprintf("maximum deviation:%.6f\n", err1);

fprintf("\nSolving the system of nonlinear equations in terms of the alternating set:\n")
fprintf("crossing point set: %.6f, %.6f, %.6f, %.6f\n", x2(1), x2(2), x2(3), x2(4));
fprintf("fitting polynomial: %.6fx^2+%.6fx+%.6f\n", c2(3), c2(2), c2(1));
fprintf("maximum deviation:%.6f\n", err2);

figure; px = linspace(-1, 1, 100);
plot(px, f(px), px, g1(px), px, g2(px));
legend("original-exp(x)", "Remez algorithm", "cross-ver");
title("fitting result");