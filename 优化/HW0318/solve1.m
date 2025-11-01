tic;
load('./data/A.mat', 'A');
load('./data/b.mat', 'b');
% b = A * u; 
x0 = zeros(1024, 1);
sig = 0.1; miu = 1;
f = @(x) lasso(A, b, sig, miu, x);
df = @(x) dflasso(A, b, sig, miu, x);
dir = @(df, x) direction(df, x);
ls = @(f, df, d, x) armijo(f, df, d, x, 6.25e-2, 0.325, 0.6);
opt_sol = gradient_descent(f, df, ls, dir, x0, 1e-6, 1000000);

fprintf('f value at Optimal solution: %.6f\n', f(opt_sol));
fprintf('norm of gradient at Optimal solution: %.6f\n', norm(df(opt_sol)));
fprintf('time spent: %f sec\n', toc);

function opt_sol = gradient_descent(f, df, ls, dir, x0, tol, maxiter)
    x1 = x0; id = 1; nd = norm(df(x1));
    f_list = zeros(1, maxiter);
    a_list = zeros(1, maxiter);
    f_list(id) = f(x1); flag = 1;
    ngl = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
    while (nd > tol && id < maxiter)
        d = dir(df, x0);
        alpha = ls(f, df, d, x0);
        x1 = x0 + alpha * d;
        x0 = x1; id = id + 1;
        f_list(id) = f(x1);
        a_list(id - 1) = alpha;
        nd = norm(df(x1));
        if (f_list(id - 1) == f_list(id))
            break;
        end
        if (nd < ngl(flag))
            fprintf("gradient level 1e%d reached \n", log10(ngl(flag)));
            fprintf('f value: %.6f norm of gradient: %.6f\n', f_list(id), nd);
            fprintf('iteration %d times, time spent: %f sec\n\n', id, toc);
            flag = flag + 1;
        end
    end
    fprintf('Number of iterations: %d\n', id);
    figure;
    plot(1:id - 1, log10(f_list(1:id - 1) - f_list(id)));
    title('f value'); ylabel('log10 f value over the lowest value');
    xlabel('iter time');
    figure;
    plot(1:id - 1, log10(a_list(1:id - 1)));
    title('log version of step size for each turn')
    ylabel("log10 of step size"); xlabel("iteration turns");
    if (maxiter == id)
        warning('maximum number of iterations reached! solution may not be optimal!');
    end
    opt_sol = x1;
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

% function [a,t] = wolfe(f, df, d, x, a0)
%     a = a0; t=0; c1=0.09; c2=0.3;
%     f0 = f(x); df0 = dot(df(x), d);
%     secant = (f(x + a * d) - f0) / a; dfx = dot(df(x + a * d), d);
%     while (secant > df0 * c1 || dfx < df0 * c2)
%         use the quadratic interpolation (f(0),f'(0),f(a)).
%         update a as the minimum point of quadratic function.
%         a = df0 * a / (2 * (df0 - secant));
%         secant = (f(x + a * d) - f0) / a; dfx = dot(df(x + a * d), d);
%         if (a < 1e-6 && secant <= 0)
%             t=1;
%             break;
%         end
%     end
% end

function d = direction(df, x)
    d = -df(x);
    % you can change the logic here to choose
    % the direction you like. For convience,
    % we just use the negative gradient.
    if (dot(d, df(x)) >= 0)
        warning('unsuccessful gradient descent direction');
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
