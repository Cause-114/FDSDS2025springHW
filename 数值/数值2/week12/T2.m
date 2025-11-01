f = @(t, y) exp(t);
actual = @(t) exp(t);
y0 = 1; t0 = 0; tk = 100;
h_list = (tk-t0)./(100:50:1000);
e_list = zeros(length(h_list), 1);
for i = 1:length(h_list)
    res = RK4(f, y0, t0, tk, h_list(i));
    h_list(i) = log(h_list(i));
    e_list(i) = log(abs(actual(tk) - res(end)));
end
[a, b] = least_squares(h_list, e_list);
fprintf('The slope of the log-log plot is %f\n', a);
figure;
plot(h_list, e_list, h_list, a * h_list + b);
xlabel('log(h)'); ylabel('log(error)');
title('Log-Log Plot of Error');
legend('Actual error', sprintf('y = %f*x + %f', a, b));

function res = RK4(f, y0, t0, tk, h)
    L = round((tk - t0) / h);
    res = zeros(L + 1, 1);
    res(1) = y0;
    for i = 2:L + 1
        k1 = h * f(t0, res(i - 1));
        k2 = h * f(t0 + h / 2, res(i - 1) + k1 / 2);
        k3 = h * f(t0 + h / 2, res(i - 1) + k2 / 2);
        k4 = h * f(t0 + h, res(i - 1) + k3);
        res(i) = res(i - 1) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        t0 = t0 + h;
    end
end

function [a, b] = least_squares(x, y)
    n = length(x);
    x_ = mean(x); y_ = mean(y);
    a = (dot(x, y) - n * x_ * y_) / (dot(x, x) - n * x_ * x_);
    b = (y_ - a * x_);
end

% function res= Euler(f, y0, t0, tk, h)
%     L = round((tk - t0) / h);
%     res = zeros(L + 1, 1);
%     res(1) = y0;
%     for i = 2:L + 1
%         res(i) = res(i - 1) + h * f(t0, res(i - 1));
%         t0 = t0 + h;
%     end
% end
