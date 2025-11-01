% In this code, I want to compute fitting g(x)=ax^2+bx+c and
% f(x)=g(x)+e^x by solving equations using Newton method.
% f1 is equal to f(1)+f(-1)=0, similarly, f2: f(x)=f(1), f5: f(x)+f(y)=0.
% f3: f'(x)=0, f4: f'(y)=0.
% where -1,x,y,1 are the crossing points set.

% Notice that for this specific problem, we assume that two boundary points are among the 
% crossing points set, and we listed all the equation of f and df. So when the function  
% is more complex or the dimension becomes lager, this method is not suitable.
function [x, c, err] = cross_ver(tol)
    e = exp(1);
    f1 = @(a, c) a + c + (e + 1 / e) / 2;
    f2 = @(a, b, x) exp(x) + a * x ^ 2 + b * x - e - a - b;
    f3 = @(a, b, x) exp(x) + 2 * a * x + b;
    f4 = @(a, b, y) exp(y) + 2 * a * y + b;
    f5 = @(a, b, c, x, y) exp(x) + exp(y) + a * (x ^ 2 + y ^ 2) + b * (x + y) + 2 * c;
    f = @(v) [f1(v(1), v(3)), f2(v(1), v(2), v(4)), f3(v(1), v(2), v(4)), ...
                  f4(v(1), v(2), v(5)), f5(v(1), v(2), v(3), v(4), v(5))];
    df = @(v)[1, v(4) ^ 2 - 1, 2 * v(4), 2 * v(5), v(4) ^ 2 + v(5) ^ 2; ...
                  0, v(4) - 1, 1, 1, v(4) + v(5); 1, 0, 0, 0, 2; ...
                  0, 2 * v(1) * v(4) + v(2) + exp(v(4)), exp(v(4)) + 2 * v(1), ...
                  0, 2 * v(1) * v(4) + v(2) + exp(v(4)); ...
                  0, 0, 0, 2 * v(1) + exp(v(5)), 2 * v(1) * v(5) + v(2) + exp(v(5))];
    v = [1, 1, 0.5, 0, 0.5];
    v = newton_solve(v, f, df, tol);
    x = [-1, v(4), v(5), 1]; c = [-v(3), -v(2), -v(1)];
    err = abs(v(1) + v(2) + v(3) + e);
end

function root = newton_solve(initial_guess, f, df, tolerance)
    iter_time = 0; fr = f(initial_guess);
    root = initial_guess;
    while (norm(fr) > tolerance)
        root = root - fr / df(root);
        fr = f(root);
        iter_time = iter_time + 1;
    end
    fprintf("total %d iterations, by newton method.\n", iter_time);
end
