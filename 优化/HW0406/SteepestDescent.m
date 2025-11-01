function x0 = SteepestDescent(tol, maxiter)
    fprintf("\n\n\nSteepest descent algorithm starts...\n");
    tic;
    load('./data/A.mat', 'A');
    load('./data/b.mat', 'b');
    x0 = zeros(1024, 1); sig = 0.5; miu = 1;
    f = @(x) lasso(A, b, sig, miu, x);
    df = @(x) dflasso(A, b, sig, miu, x);
    ls = @(f, df, d, x) armijo(f, df, d, x, 6.25e-2, 0.325, 0.6);
    id = 1; nd = norm(df(x0));
    flag = 1; ngl = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
    while (nd > tol && id < maxiter)
        d = -df(x0);
        alpha = ls(f, df, d, x0);
        x0 = x0 + alpha * d;
        id = id + 1; nd = norm(d);
        if (nd < ngl(flag))
            fprintf("gradient level 1e%d reached \n", log10(ngl(flag)));
            fprintf("iteration %d times\n", id);
            flag = flag + 1;
        end
    end
    fprintf("time spent: %f sec\n", toc);
    fprintf("Number of iterations: %d\n", id);
    if (maxiter == id)
        warning('maximum number of iterations reached! solution may not be optimal!');
    end
    fprintf("f value at Optimal solution: %.6f\n", f(x0));
    fprintf("norm of gradient at Optimal solution: %.6f\n", nd);
    fprintf("Steepest descent algorithm ends...\n\n\n");
end

function a0 = armijo(f, df, d, x, a0, c1, beta)
    f0 = f(x); df0 = dot(df(x), d);
    fx = f(x + a0 * d);
    while (fx > f0 + c1 * a0 * df0)
        a0 = a0 * beta; fx = f(x + a0 * d);
        if (a0 < 1e-10 && fx <= f0)
            break;
        end
    end
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
