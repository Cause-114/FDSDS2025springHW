tic;
load('./data/A.mat', 'A');
load('./data/b.mat', 'b');
sig = 0.1; miu = 0.01; a0 = 1e-3;
f = @(x) lasso(A, b, sig, miu, x);
df = @(x) dflasso(A, b, sig, miu, x);
ls = @(f, df, d, x) wolfe(f, df, d, x, a0);
dir = @(df, x) direction(df, x);
x0 = zeros(size(A, 2), 1);
opt_sol = gradient_descent(f, df, ls, dir, x0, 1e-6, 1000000);
fprintf('f value at Optimal solution: %.6f\n', f(opt_sol));
fprintf('norm of gradient at Optimal solution: %.6f\n', norm(df(opt_sol)));
fprintf('time spent: %f sec\n',toc);
function opt_sol = gradient_descent(f, df, ls, dir, x0, tol, maxiter)
    x1 = x0; id = 1;
    f_list = zeros(1, maxiter);
    f_list(id) = f(x1);
    a_list= zeros(1,maxiter);
    while (norm(df(x1)) > tol && id < maxiter)
        d = dir(df, x0); alpha = ls(f, df, d, x0);
        x1 = x0 + alpha * d;
        x0 = x1; id = id + 1;
        f_list(id) = f(x1);
        a_list(id-1)=alpha;
        if(f_list(id-1)==f_list(id))
            break;
        end
    end
    fprintf('Number of iterations: %d\n', id); figure;
    plot(1:id-1, log10(f_list(1:id-1)-f_list(id))); 
    title('f value');ylabel('log10 f value over the lowest value');
    xlabel('iter time');figure;
    plot(1:id-1,log10(a_list(1:id-1)))
    if (maxiter == id)
        warning('maximum number of iterations reached! solution may not be optimal!');
    end
    opt_sol = x1;
end

% function a = armijo(f, df, d, x, a0)
%     a = a0; f0 = f(x); df0 = dot(df(x), d);
%     c1 = 0.085; beta = 0.35; fx = f(x + a * d);
%     while (fx > f0 + c1 * a * df0)
%         a = a * beta; fx = f(x + a * d);
%         if (a < 1e-10 && fx <= f0)
%             break;
%         end
%     end
% end

function a = wolfe(f, df, d, x, a0)
    a = a0; c1 = 0.09; c2 = 0.3;
    f0 = f(x); df0 = dot(df(x), d);
    secant = (f(x + a * d) - f0) / a; dfx = dot(df(x + a * d), d);
    while (secant > df0 * c1 || dfx < df0 * c2)
        % use the quadratic interpolation (f(0),f'(0),f(a)).
        % update a as the minimum point of quadratic function.
        a = df0 * a / (2 * (df0 - secant));
        secant = (f(x + a * d) - f0) / a; dfx = dot(df(x + a * d), d);
        if (a < 1e-6 && secant <= 0)
            break;
        end
    end
end

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
