% x stands for -pH, y stands for -pOH, z stands for -pHCO3-, and w stands for -pCO3-2.
% by transforming the equations into a more linear form, we can simplify the caculating process.
f1 = @(x, y) x + y + 14;
f2 = @(x, z) x + z + 13.76 - log10(375);
f3 = @(x, z, w) x + w - z + 10.3;
f4 = @(x, y, z, w) 10 ^ x - 10 ^ y - 10 ^ z - 2 * 10 ^ w;
% In the process of newton iteration, in order to use the same form of x_{k+1}=x_k-f(x_k)/df(x_k),
% we transpose the vector f and matrix df.
f = @(vari) [f1(vari(1), vari(2)), f2(vari(1), vari(3)), ...
                 f3(vari(1), vari(3), vari(4)), f4(vari(1), vari(2), vari(3), vari(4))];
df = @(vari)[1, 1, 1, 10 ^ vari(1) * log(10); ...
                 1, 0, 0, -10 ^ vari(2) * log(10); ...
                 0, 1, -1, -10 ^ vari(3) * log(10); ...
                 0, 0, 1, -2 * 10 ^ vari(4) * log(10)];
fprintf("\tpH\tpOH\tpHCO3-\tpCO3-2\n");
initial_guess = [-6, -8, -6, -10]; tolerance = 10 ^ -12;
root = newton_solve(initial_guess, f, df, tolerance);
disp(-root);
root1 = broyden_solveg(initial_guess, f, tolerance, root);
disp(-root1);
root2 = broyden_solveb(initial_guess, f, tolerance, root);
disp(-root2);

function r = newton_solve(initial_guess, f, df, tolerance)
    iter_time = 1; r = initial_guess;
    dr = f(r) / df(r);
    while norm(dr) > tolerance
        r = r - dr;
        dr = f(r) / df(r);
        iter_time = iter_time + 1;
    end
    fprintf('Number of iterations: %d, by newton method.\n', iter_time);
end

% function root = broyden_solve1(initial_guess, f, df, tolerance)
%     % this is the 'good' Broyden method (by meshing the sherman-morrison formula).
%     J = df(initial_guess);
%     x0 = initial_guess; fx0 = f(x0);
%     x1 = x0 - fx0 / J; fx1 = f(x1);
%     iter_time = 1; J = inv(J);
%     while norm(x1 - x0) > tolerance
%         iter_time = iter_time + 1;
%         dx = x1 - x0; df = fx1 - fx0;
%         J = J + J * dx' * (dx - df * J) / (df * J * dx');
%         x0 = x1; fx0 = fx1;
%         x1 = x0 - fx0 * J; fx1 = f(x1);
%     end
%     fprintf('Number of iterations: %d, by good broyden.\n', iter_time);
%     root = x1;
% end
function x0 = broyden_solveg(initial_guess, f, tolerance, root)
    % In this version, we don't use the Jacobian matrix.
    % this is the 'good' Broyden method (by meshing the sherman-morrison formula).
    J = eye(length(initial_guess));
    x0 = initial_guess; fx0 = f(x0);
    iter_time = 1; dx =- fx0 * J;
    conver= zeros(1,50);conver(1)=norm(x0-root);
    while norm(dx) > tolerance
        x0 = x0 + dx; df = f(x0) - fx0;
        J = J + J * dx' * (dx - df * J) / (df * J * dx');
        fx0 = fx0 + df; dx =- fx0 * J;
        iter_time = iter_time + 1;
        conver(iter_time)=norm(x0-root);
    end
    fprintf('Number of iterations: %d, by good broyden.\n', iter_time);
    figure;
    plot(1:iter_time,log10(conver(1:iter_time)));
    title('Convergence of Broyden method(Good)');
end

% function root = broyden_solve2(initial_guess, f, df, tolerance)
%     % this is the 'bad' Broyden method (by using the formula to the inverse of the Jacobian).
%     J = df(initial_guess);
%     x0 = initial_guess; fx0 = f(x0);
%     x1 = x0 - fx0 / J; fx1 = f(x1);
%     iter_time = 1; J = inv(J);
%     while norm(x1 - x0) > tolerance
%         iter_time = iter_time + 1; df = fx1 - fx0;
%         J = J + df' * (- fx1 * J) / dot(df, df);
%         x0 = x1; fx0 = fx1;
%         x1 = x0 - fx0 * J; fx1 = f(x1);
%     end
%     fprintf('Number of iterations: %d, by bad broyden.\n', iter_time);
%     root = x1;
% end

function x0 = broyden_solveb(initial_guess, f, tolerance, root)
    % In this version, we don't use the Jacobian matrix.
    % this is the 'bad' Broyden method (by using the formula to the inverse of the Jacobian).
    J = eye(length(initial_guess));
    x0 = initial_guess; fx0 = f(x0);
    iter_time = 1;
    conver= zeros(1,50);conver(1)=norm(x0-root);
    while (norm(fx0 * J) > tolerance)
        x0 = x0 - fx0 * J; df = f(x0) - fx0;
        J = J + df' * (- (fx0 + df) * J) / dot(df, df);
        fx0 = fx0 + df;
        iter_time = iter_time + 1;
        conver(iter_time)=norm(x0-root);
    end
    fprintf('Number of iterations: %d, by bad broyden.\n', iter_time);
    figure;
    plot(1:iter_time,log10(conver(1:iter_time)));
    title('Convergence of Broyden method(Bad)');
end
% f=@(x) x^3-(10^-14+375*10^-13.76)*x-375*2*10^-24.06;
% df=@(x) 3*x^2-(10^-14+375*10^-13.76);
% root=newton_solve(1,f,df,10^-15);
% fprintf('Root of the function is: log10(root) = %f\n',log10(root));
