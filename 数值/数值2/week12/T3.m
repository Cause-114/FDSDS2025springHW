f = @(x) derivative(x(1), x(2), x(3));
t0 = 0; tk = 100; h = 0.01;
y0 = [1, 0, 100];
res = RK4(f, y0, t0, tk, h);
figure;
plot(t0:h:tk, res(:, 1), t0:h:tk, res(:, 2), t0:h:tk, res(:, 3));
title("Concentrations over time"); xlabel("t");
ylabel("concentrations");
legend("bacteria", "detritus", "substrate")
function res = RK4(f, y0, t0, tk, h)
    L = round((tk - t0) / h);
    res = zeros(L + 1, 3);
    res(1, :) = y0;
    for i = 2:L + 1
        k1 = h * f(res(i - 1, :));
        k2 = h * f(res(i - 1, :) + k1 / 2);
        k3 = h * f(res(i - 1, :) + k2 / 2);
        k4 = h * f(res(i - 1, :) + k3);
        res(i, :) = res(i - 1, :) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
end

function res = derivative(X, C, S)
    miu_m = 10; K = 10; Ks = 10;
    kd = 0.1; ke = 0.1; kh = 0.1;
    res = zeros(1, 3);
    res(1) = miu_m * (1 - X / K) * (S / (Ks + S)) * X - kd * X - ke * X;
    res(2) = kd * X - kh * C;
    res(3)=-res(1)-res(2);
end
