% We have known that solving a linear system costs O(n^3) time, 
% so if we iterates many time, it's not practical for large n.
% Here we use interplotation method to solve the linear system,
% which instead costs O(n^2) time. However, the expression form
% of interplotation results is not easy to find 1st\2nd derivatives,
% which makes is difficult to use newton method to update the new crossing
% points. Therefore, we use bisection method to update the crossing points.
% And it can be proved that given specific tolerances, bisection method
% will stop in limited number of iterations. So, each nodes cost O(n) time
% to caculate derivatives, which makes the total time complexity O(n^2).
% As shown above each turn we costs O(n^2).

% What calls for special attention is the output form of polynomials,
% it looks like: p(x) = c1(1) + c1(2)*(x-x(1)) + c1(3)*(x-x(1))*(x-x(2)) +...
function [x, c1,err] = remez(f, df, x, a, b, tol)
    n = length(x);iter_time=1;
    while true
        % disp(x);
        [c1, E] = inter(f, x);
        nx = update_x(f, df, x, c1, a, b, tol);
        t = max(abs(arrayfun(@(i) fp(f, x, c1, nx(i)), 1:n)));
        if (t - abs(E) < tol)
            err = t;
            break;
        end
        x = nx;iter_time=iter_time+1;
    end
    fprintf("total %d iterations, by remez algorithm.\n",iter_time);
end

% caculate the derivative of f(v)-p(v) at v
% p(v) is the form displayed above.
function dfv = dfp(df, x, c1, v)
    dfv = c1(2) - df(v); n = length(c1);
    term1 = 1; term2 = v - x(1);
    for i = 3:n
        term1 = term1 * (v - x(i - 1)) + term2;
        dfv = dfv + term1 * c1(i);
        term2 = term2 * (v - x(i - 1));
    end
end

% caculate the value of f(v)-p(v) at v
% p(v) is the form displayed above.
function fv = fp(f, x, c1, v)
    fv = c1(1) - f(v); pd = 1;
    n=length(x);
    for j = 2:n - 1
        pd = pd * (v - x(j - 1));
        fv = fv + c1(j) * pd;
    end
end

function m = bisection(f, l, r, tol)
    sgn = sign(f(l));
    while (r - l > tol)
        m = (l + r) / 2; fm = f(m);
        if abs(fm) < tol
            break;
        elseif fm * sgn > 0
            l = m;
        else
            r = m;
        end
    end
end

% When given crossing points x(1), x(2),..., x(n), we find the coefficients of the 
% polynomial p(x) and the deviation E = f(x(1)) - p(x(1)) by using the interplotation method
function [c1, E] = inter(f, x)
    n = length(x);
    c1 = f(x(1:n - 1)); c2 = mod(1:n - 1, 2) * 2 - 1;
    pd = 1; p1 = c1(1); p2 = c2(1);
    for i = 2:n - 1
        tmp = (c1(i) - c1(i - 1)) / (x(i) - x(1));
        tmp2 = (c2(i) - c2(i - 1)) / (x(i) - x(1));
        for j = i + 1:n - 1
            tmp1 = (c1(j) - c1(j - 1)) / (x(j) - x(j - i + 1));
            tmp3 = (c2(j) - c2(j - 1)) / (x(j) - x(j - i + 1));
            c1(j - 1) = tmp; tmp = tmp1;
            c2(j - 1) = tmp2; tmp2 = tmp3;
        end
        c1(n - 1) = tmp; c2(n - 1) = tmp2;
        pd = pd * (x(n) - x(i - 1));
        p1 = p1 + c1(i) * pd;
        p2 = p2 + c2(i) * pd;
    end
    E = (p1 - f(x(n))) / (p2 + 1);
    c1 = c1 - E * c2;
end

% First find the roots of f(v)-p(v) between x(i) and x(i+1)
% Then in the interval of each roots, there must lies a local min(max) points.
% For the boundary points, we specially handle them.
function nx = update_x(f, df, x, c1, a, b, tol)
    n = length(x);
    g = @(v) fp(f, x, c1, v);
    rt = arrayfun(@(i) bisection(g, x(i), x(i + 1), tol), 1:n - 1);
    nx = zeros(1, n); g = @(v) dfp(df, x, c1, v);
    nx(2:n - 1) = arrayfun(@(i) bisection(g, rt(i - 1), rt(i), tol), 2:n - 1);
    if (g(a) * g(rt(1)) > 0)
        nx(1) = a;
    else
        nx(1) = bisection(g, a, rt(1), tol);
    end
    if (g(b) * g(rt(n - 1)) > 0)
        nx(n) = b;
    else
        nx(n) = bisection(g, a, rt(n - 1), tol);
    end
end
