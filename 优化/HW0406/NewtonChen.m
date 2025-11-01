% I recommend you not to set tolerence to a quite small value
% For my test, it can only reach the 1e-11.75 level of gradient norm.
% Here is the recommended value: tol=1e-10, maxiter=10.
% You can adjust the value according to your need.
function x0 = NewtonChen(tol, maxiter)
    fprintf("\n\n\nChen's code for Newton Method starts...\n")
    tic;
    load('./data/A.mat', 'A');
    load('./data/b.mat', 'b');
    x0 = zeros(1024, 1);
    sig = 0.5; miu = 1;
    f = @(x) lasso(A, b, sig, miu, x);
    df = @(x) dflasso(A, b, sig, miu, x);
    dir = @(df, x) direction(df, x, A, sig, miu);
    ls = @(f, df, d, x) armijo(f, df, d, x, 1, 0.1, 0.6);
    id = 1; nd = norm(df(x0));
    while (nd > tol && id < maxiter)
        d = dir(df, x0);
        alpha = ls(f, df, d, x0);
        x0 = x0 + alpha * d;
        id = id + 1;
        nd = norm(df(x0));
    end
    fprintf("time spent: %f sec\n", toc);
    fprintf("Number of iterations: %d\n", id);
    if (maxiter == id)
        warning("maximum number of iterations reached! solution may not be optimal!");
    end
    fprintf("f value at Optimal solution: %.6f\n", f(x0));
    fprintf("log10-version of norm of gradient at Optimal solution:%.2f", log10(nd));
    fprintf("\nChen's code for Newton Method ends...\n\n\n\n");
end

function a0 = armijo(f, df, d, x, a0, c1, beta)
    f0 = f(x); df0 = dot(df(x), d);
    fx = f(x + a0 * d);
    while (fx > f0 + c1 * a0 * df0)
        a0 = a0 * beta; fx = f(x + a0 * d);
    end
end

function d = direction(df, x, A, sig, miu)
    H = A' * A + miu * diag((abs(x) <= sig) ./ sig);
    H = H + eye(size(H)) * 0;
    d = -H \ df(x);
end

function res = lasso(A, b, sig, miu, x)
    abs_x = abs(x);
    quadratic_part = (abs_x <= sig) .* (x .^ 2 / (2 * sig));
    linear_part = (abs_x > sig) .* (abs_x - sig / 2);
    penalty = sum(quadratic_part + linear_part) * miu;
    res = norm(A * x - b) ^ 2/2 + penalty;
end

function dfres = dflasso(A, b, sig, miu, x)
    abs_x = abs(x);
    gradient_penalty = (abs_x > sig) .* sign(x) + (abs_x <= sig) .* (x / sig);
    dfres = gradient_penalty * miu + A' * (A * x - b);
end
