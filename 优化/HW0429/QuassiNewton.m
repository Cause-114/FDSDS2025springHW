main();

function main()
    if ~exist('./a1a.txt', 'file')
        error('File a1a.txt not found. Please check the file path.');
    end
    [b, A] = libsvmread('./a1a.txt');
    miu = 1e-3; % initialise parameters A,b,miu
    f = @(x) logisticRegression(miu, A, b, x);
    df = @(x) grad_LR(miu, A, b, x); % initialize function handle of f,df,lineSearch.
    ls = @(d, x) lineSearch(f, df, d, x, 1.5, 0.01, 0.99);
    n = size(A, 2); x0 = zeros(n, 1);
    maxIter = 5000; tol = 1e-8; % set initial guess, tolerance and max iteration.
    % The session below aims to compare the results get by DFP and BFGS, if you only want
    % to justify one of my self-implemented algorithm, you can delete the other part,
    
    % if you want to make pictures, remember to add f in the fucntion handle.
    % and use the code that are commented at the bottom of the file.
    x1 = BFGS(df, ls, x0, maxIter, tol);
    x2 = DFP(df, ls, x0, maxIter, tol);
    fprintf("\nNorm of optimal point difference between BFGS and DFP: %.8g\n", norm(x1 - x2));
    fprintf("Gradient norm: %g(BFGS), %g(DFP)\n", norm(df(x1)), norm(df(x2)));
    fprintf("F-value at optimal solution: %.8g(BFGS), %.8g(DFP)\n\n", f(x1), f(x2));
    % save("./x1.mat", "x1");
    % SD(f, df, ls, x0, maxIter, tol);
end

% logistic regression function:
% 1/m*Sigma_{i=1}^m log(1+exp(-bi.*ai*x)) + miu*|x|^2
function val = logisticRegression(miu, A, b, x)
    m = length(b);
    val = sum(log(1 + exp(-b .* A * x)));
    val = val / m + miu * dot(x, x);
end

% gradient of logistic regression function:
% 1/m*Sigma_{i=1}^m ((exp(-bi.*ai*x)/(1+exp(-bi.*ai*x))) * -bi .* ai) + 2*miu*x
function g = grad_LR(miu, A, b, x)
    m = length(b);
    g = sum((1 ./ (1 + exp(-b .* A * x)) - 1) .* b .* A);
    g = g' / m + 2 * miu * x;
end

function a = lineSearch(f, df, d, x, a, c, beta)
    f0 = f(x); df0 = dot(df(x), d); minStep = 1e-8;
    while (a > minStep && (f(x + a * d) > f0 + c * a * df0||...
            df0 >= dot(df(x + a * d), d)))
    % while make sure f value descent, we further require df_{k+1}^Td>df_k^Td
        a = a * beta;
    end
    if (a <= minStep)
        warning("Line search failed. Using minimum step size instead.");
    end
end

function x0 = BFGS(df, ls, x0, maxIter, tol)
    y0 = df(x0); H = eye(length(x0)); id = 1;
    while (norm(y0) > tol && maxIter > id)
        s =- H * y0; d = ls(s, x0);
        s = s * d; x0 = x0 + s;
        y = df(x0) - y0;
        % note that (I-\rho_k s_k y_k^T) H_k(I-\rho_k y_k s_k^T)+\rho_k s_k s_k^T
        % so the update terms contains Hys'(and its transpose) and ss'
        % In this way, we avoid the O(n^3) cost for matrix mutiplication.
        t2 = s' * y; t1 = H * y * (s' / t2);
        H = H - (t1 + t1') + s * (s' * ((t2 + y' * H * y) / (t2 ^ 2)));
        y0 = y0 + y; id = id + 1;
    end
    if (id == maxIter)
        fprintf("Warning of reach max iteration turns at BFGS.\n");
    end
    fprintf("Total %d iterations to get %g accuracy by my self-implemented BFGS\n", id, tol);
end

function x0 = DFP(df, ls, x0, maxIter, tol)
    y0 = df(x0); H = eye(length(x0)); id = 1;
    while (norm(y0) > tol && maxIter > id)
        s =- H * y0; d = ls(s, x0);
        s = s * d; x0 = x0 + s;
        y = df(x0) - y0;
        % nothing special, just the normal upate formula
        % but remember to avoid the O(n^3) cost for matrix mutiplication.
        t1 = H * y;
        H = H - t1 * (t1' / (y' * H * y)) + s * (s' / (s' * y));
        y0 = y0 + y; id = id + 1;
    end
    if (id == maxIter)
        fprintf("Warning of reach max iteration turns at DFP.\n");
    end
    fprintf("Total %d iterations to get %g accuracy by my self-implemented DFP\n", id, tol);
end

%% Bellow part is just for piture make.
% It asks for additional costs of storage spaces and need f handles passed
% to draw the relationship between f and iteration turns. So when you want
% to use it, you need to modify interface at line 16-17 to add f handle.

% function x0 = BFGS(f, df, ls, x0, maxIter, tol)
%     y0 = df(x0); H = eye(length(x0)); id = 1; tol = log10(tol);
%     norm_y0 = zeros(maxIter, 1); norm_y0(1) = log10(norm(y0));
%     y_list = zeros(maxIter, 1); y_list(1) = f(x0);
%     while (norm_y0(id) > tol && maxIter > id)
%         s =- H * y0; d = ls(s, x0);
%         s = s * d; x0 = x0 + s;
%         y = df(x0) - y0;
%         t2 = s' * y; t1 = H * y * (s' / t2);
%         H = H - (t1 + t1') + s * (s' * ((t2 + y' * H * y) / (t2 ^ 2)));
%         y0 = y0 + y; id = id + 1;
%         norm_y0(id) = log10(norm(y0));
%         y_list(id) = f(x0);
%     end
%     if (id == maxIter)
%         fprintf("Warning of reach max iteration turns at BFGS.\n");
%     end
%     fprintf("Total %d iterations to get %g accuracy by my self-implemented BFGS\n", id, 10^tol);
%     figure; plot(1:id, norm_y0(1:id)); title("Norm of gradient at each iteration(BFGS)");
%     xlabel("Iteration"); ylabel("log10(Norm of gradient)");
%     figure; plot(1:id-1, log10(y_list(1:id-1)-y_list(id)));
%     title("log-version of F-value over optimal solution at each iteration(BFGS)");
%     xlabel("Iteration"); ylabel("log10(F-value difference)");
% end

% function x0 = DFP(f, df, ls, x0, maxIter, tol)
%     y0 = df(x0); H = eye(length(x0)); id = 1; tol = log10(tol);
%     norm_y0 = zeros(maxIter, 1); norm_y0(1) = log10(norm(y0));
%     y_list = zeros(maxIter, 1); y_list(1) = f(x0);
%     while (norm_y0(id) > tol && maxIter > id)
%         s =- H * y0; d = ls(s, x0);
%         s = s * d; x0 = x0 + s;
%         y = df(x0) - y0;
%         t1 = H * y;
%         H = H - t1 * (t1' / (y' * H * y)) + s * (s' / (s' * y));
%         y0 = y0 + y; id = id + 1;
%         norm_y0(id) = log10(norm(y0));
%         y_list(id) = f(x0);
%     end
%     if (id == maxIter)
%         fprintf("Warning of reach max iteration turns at DFP.\n");
%     end
%     fprintf("Total %d iterations to get %g accuracy by my self-implemented DFP\n", id, 10^tol);
%     figure; plot(1:id, norm_y0(1:id)); title("Norm of gradient at each iteration(DFP)");
%     xlabel("Iteration"); ylabel("log10(Norm of gradient)");
%     figure; plot(1:id-1, log10(y_list(1:id-1)-y_list(id)));
%     title("log-version of F-value over optimal solution at each iteration(DFP)");
%     xlabel("Iteration"); ylabel("log10(F-value difference)");
% end

% function x0 = SD(f, df, ls, x0, maxiter, tol)
%     % ls = @(f, df, d, x) armijo(f, df, d, x, 1.5, 0.01, 0.99);
%     nd=zeros(maxiter,1); nf=zeros(maxiter,1);
%     id = 1; nd(id) = norm(df(x0)); nf(id) = f(x0);
%     while (nd(id) > tol && id < maxiter)
%         d = -df(x0);
%         alpha = ls(d, x0);
%         x0 = x0 + alpha * d;
%         id = id + 1; nd(id) = norm(d);
%         nf(id) = f(x0);
%     end
%     if (maxiter == id)
%         warning('maximum number of iterations reached! solution may not be optimal!');
%     end
%     fprintf("with accuracy %g, Number of iterations: %d (SD)\n", tol, id);
%     fprintf("f value at Optimal solution: %.8g (SD)\n", nf(id));
%     fprintf("norm of gradient at Optimal solution: %g (SD)\n", nd(id));
%     figure; plot(1:id, log10(nd(1:id))); title("Norm of gradient at each iteration(SD)");
%     xlabel("Iteration"); ylabel("log10(Norm of gradient)");
%     figure; plot(1:id-1, log10(nf(1:id-1)-nf(id)));
%     title("log-version of F-value over optimal solution at each iteration(SD)");
%     xlabel("Iteration"); ylabel("log10(F-value difference)");
% end
