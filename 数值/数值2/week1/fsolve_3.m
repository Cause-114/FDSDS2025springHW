% 定义方程
fun = @(x) (2 * x / (1 + x ^ 2)) - atan(x);

% % 提供一个初始猜测值。对于这个方程，可以从 x=1 开始尝试。
% x0 = 1;
% options = optimset('TolX', 1e-11, 'TolFun', 1e-11);
% % 使用 fzero 求解
% solution = fzero(fun, x0,options);

% % 显示结果
% fprintf("%.11f",solution);
x1 = 1; x2 = 3;
while (x2 - x1 > 1e-12)
    mid = (x1 + x2) / 2;
    if (fun(mid) > 0)
        x1 = mid;
    else
        x2 = mid;
    end
end

fprintf("The ans is %.12f", (x1 + x2) / 2);
