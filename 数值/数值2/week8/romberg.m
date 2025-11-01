f = @(x) exp(x);
g = @(x) x .^ 1.5;
a = 0; b = 1;
format_T(Romberg(f, a, b, 5, exp(1)-1));
format_T(Romberg(g, a, b, 5, 0.4));

function T = Romberg(f, a, b, n, res)
    T = zeros(1, n + 1); h = b - a;
    T(1) = h * (f(a) + f(b)) / 2;
    for i = 1:n
        h = h / 2; d = 4;
        tmp = T(1) / 2 + h * sum(f(a + h * (1:2:2 ^ i)));
        for j = 2:i + 1
            tmp2 = (tmp * d - T(j - 1)) / (d - 1);
            T(j - 1) = tmp; tmp = tmp2; d = d * 4;
        end
        T(i + 1) = tmp;
        e = T(1:i + 1) - res;
        format_e(e);
    end
end

function format_e(e)
    fprintf("value of error in the %d-th row:\n", length(e))
    for i = 1:length(e)
        fprintf("%.4e\t", e(i));
    end
    fprintf("\n");
end

function format_T(T)
    fprintf("\nvalue of T after %d step:\n", length(T) - 1);
    for i = 1:length(T)
        fprintf("%.4e\t", T(i));
    end
    fprintf("\n\n");
end
