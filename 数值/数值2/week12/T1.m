f = @(x) arrayfun(@(x) discontinuous_function(x), x);
filter1 = @(x) arrayfun(@(x) Gaussian(x), x);
filter2 = @(x) arrayfun(@(x) Friedrichs(x), x);

px = linspace(-10, 10, 1000);
plot(px, f(px));
hold on;
plot(px, convolve(f, filter1, px));
hold on;
plot(px, convolve(f, filter2, px));
title('Filter Comparison');
xlabel('x'); ylabel('y');
legend('f(x)', 'Gaussian Filter', 'Friedrichs mollifiers Filter');

function py = convolve(f, filter, px)
    py = zeros(size(px));
    x = linspace(-10, 10, 100);
    tau = x(2) - x(1);
    fx = f(x);
    for i = 1:length(px)
        tmp = fx .* filter(px(i) - x);
        py(i) = trapz(tau, tmp);
    end
end

function res = Friedrichs(x)
    if (abs(x) >= 1)
        res = 0;
    else
        res = 2.2523 * exp(1 / (x ^ 2 - 1));
    end
end

function res = Gaussian(x)
    res = 1 / sqrt(2 * pi) * exp(-x ^ 2/2);
end

function res = discontinuous_function(x)
    if (x < -10)
        res = 0;
    elseif (x <- 5)
        res = -sin(x);
    elseif (x < 0)
        res = cos(x);
    elseif (x < 5)
        res = asin(mod(x, 2) - 1);
    elseif (x < 10)
        res = acos(mod(x, 2) - 1);
    else
        res = 0;
    end
end
