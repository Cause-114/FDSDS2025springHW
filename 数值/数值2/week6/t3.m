x = [-4.000; -3.500; -3.000; -2.500; -2.000; -1.500; -1.000; -0.500; ...
         0.000; 0.500; 1.000; 1.500; 2.000; 2.500; 3.000; 3.500; 4.000];
y = [0.00001; 0.00726; 0.25811; 1.87629; 1.55654; 0.17209; 0.00899; 0.05511; ...
         0.24564; 0.60455; 0.89370; 1.03315; 0.51633; 0.18032; 0.04287; 0.00360; 0.00045];
f = @(c, x) c(1) * exp(-c(2) * (x - c(3)) .^ 2) + c(4) * exp(-c(5) * (x - c(6)) .^ 2);
c0 = [1.8; 1; -2.5; 1.5; 1; 1];
c = Guass_Newton(x, y, c0, f);
xi = (linspace(-4, 4, 1000))';
yi = f(c, xi);
figure;
plot(x, y, 'o', xi, yi, '-'); title("fitting result of samples");
xlabel('x'); ylabel('y'); legend('sample points', 'output fit function');
fprintf("a_1=%.6f\tb_1=%.6f\tc_1=%.6f", c(1), sqrt(c(2)), c(3));
fprintf("\na_2=%.6f\tb_2=%.6f\tc_2=%.6f\n", c(4), sqrt(c(5)), c(6));

function c = Guass_Newton(x, y, c, f)
    df1 = @(a, b, c, x) exp(-b * (x - c) .^ 2);
    df2 = @(a, b, c, x) -a * (x - c) .^ 2 .* exp(-b * (x - c) .^ 2);
    df3 = @(a, b, c, x) 2 * a * b * (x - c) .* exp(-b * (x - c) .^ 2);
    iter_time = 0;
    while (true)
        J = [df1(c(1), c(2), c(3), x), df2(c(1), c(2), c(3), x), df3(c(1), c(2), c(3), x), ...
                 df1(c(4), c(5), c(6), x), df2(c(4), c(5), c(6), x), df3(c(4), c(5), c(6), x)];
        b = f(c, x) - y; c = c - J \ b;
        if(norm(b' * J, 1) < 1e-10)
            break;    
        end
        iter_time = iter_time + 1;
    end
    fprintf("%d times of iteration\n", iter_time);
end
