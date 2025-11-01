y0 = [0; 1];
t0 = 0; tk = 13; h = 0.1;
method = {"Explicit Euler", "Implicit Euler","Trapezoidal", "RK4"};
figure_maker(method, y0, t0, tk, h);

function figure_maker(method, y0, t0, tk, h)
    figure;
    legends = cell(1, length(method) + 1);
    plot(sin(linspace(0, 2*pi, 100)), cos(linspace(0,2*pi,100))); 
    hold on; legends{1} = "exact";
    for s = 1:length(method)
        if strcmp(method{s}, "Explicit Euler")
            res = explicit_Euler(@derivative, y0, t0, tk, h);
        elseif strcmp(method{s}, "Implicit Euler")
            res = implicit_Euler(y0, t0, tk, h);
        elseif strcmp(method{s}, "Trapezoidal")
            res = trapezoidal(y0, t0, tk, h);
        elseif strcmp(method{s}, "RK4")
            res = RK4(@derivative, y0, t0, tk, h);
        end
        plot(res(1, :), res(2, :)); hold on;
        legends{s + 1} = method{s};
    end
    title(sprintf("set end time to be %.1f, step size to be %.1f", tk, h));
    legend(legends); xlabel("x"); ylabel("x'");
end

function res = explicit_Euler(f, y0, t0, tk, h)
    L = round((tk - t0) / h);
    res = zeros(length(y0), L + 1);
    res(:, 1) = y0;
    for i = 2:L + 1
        res(:, i) = res(:, i - 1) + h * f(res(:, i - 1));
    end
end

function res = implicit_Euler(y0, t0, tk, h)
% Note that if we use a complex derivative function, the process of finding roots will 
% be complicated for the implict method, so here I just merge derivative into the function.
    B = [1, -h; h, 1] \ [1, 0; 0, 1];
    L = round((tk - t0) / h);
    res = zeros(length(y0), L + 1);
    res(:, 1) = y0;
    for i = 2:L + 1
        res(:, i) = B * res(:, i - 1);
    end
end

function res = trapezoidal(y0, t0, tk, h)
% Note that if we use a complex derivative function, the process of finding roots will 
% be complicated for the implict method, so here I just merge derivative into the function.
    B = [1, -h / 2; h / 2, 1] \ [1, h / 2; -h / 2, 1];
    L = round((tk - t0) / h);
    res = zeros(length(y0), L + 1);
    res(:, 1) = y0;
    for i = 2:L + 1
        res(:, i) = B * res(:, i - 1);
    end
end

function res = RK4(f, y0, t0, tk, h)
    L = round((tk - t0) / h);
    res = zeros(length(y0), L + 1);
    res(:, 1) = y0;
    for i = 2:L + 1
        k1 = h * f(res(:, i - 1));
        k2 = h * f(res(:, i - 1) + k1 / 2);
        k3 = h * f(res(:, i - 1) + k2 / 2);
        k4 = h * f(res(:, i - 1) + k3);
        res(:, i) = res(:, i - 1) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
end

function df = derivative(x)
    k = 1; m = 1;
    x = x(:);
    df = [0, 1; -k / m, 0] * x;
end
