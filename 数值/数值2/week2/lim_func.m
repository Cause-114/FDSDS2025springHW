f = @(x)(x ^ 3 * exp(2 * x) - 1 - 3 * log(x)) / x;
df = @(x) exp(2 * x) * (2 * x ^ 4 + 2 * x ^ 3) - 2 + 3 * log(x);

l = 0.5; r = 1;
while (r - l > 1e-10)
    m = (l + r) / 2;
    if (abs(df(m)) < 1e-10)
        break;
    elseif (df(m) > 0)
        r = m;
    else
        l = m;
    end
end
m = (l + r) / 2;
fprintf('The root is: %.10f,\nand the inf of f(x) at the root is: %f\n', m, f(m));
