f = @(x) x .^ 3 * exp(x);
x0 = 1;
format_e(Richardson_derivative(f, x0, 7),4*exp(1));

function T = Richardson_derivative(f, x0, n)
    T = zeros(1, n + 1); h = 0.5;
    T(1) = (f(x0 + h) - f(x0 + h)) / (2 * h);
    for i = 1:n
        h = h / 2; d = 4;
        tmp = (f(x0 + h) - f(x0 - h)) / (2 * h);
        for j = 2:i + 1
            tmp2 = (tmp * d - T(j - 1)) / (d - 1);
            T(j - 1) = tmp; tmp = tmp2; d = d * 4;
        end
        T(i + 1) = tmp;
    end
end

function format_e(T,res)
    T=T-res; n=length(T)-1;
    figure;
    plot(1:n+1,log10(abs(T)));
    xlabel(sprintf("the order of extrapolation in the %d-row",n));
    ylabel("$log_{10}(|T-res|)$","Interpreter","latex");
    title(sprintf("Richardson extrapolation error after %d steps",n));
end
