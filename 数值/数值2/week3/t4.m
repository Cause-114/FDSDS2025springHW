x = linspace(0, 3 * pi, 100);
% the sample for plotting curve.
for inter_number = 3:2:9
    % bellow is the sample for interpolation.
    inter_point = linspace(0, 2 * pi, inter_number);
    f_list = sin(inter_point);
    % res stands for the interpolated function.
    res = interplot_Newton(f_list, inter_point);
    figure
    plot(x, arrayfun(res, x), "-")
    hold on
    plot(inter_point, f_list, "o")
    hold on
    plot(x, sin(x), "--")
    legend('Interpolated Curve', 'Data Points', 'original curve')
    title(sprintf("number of interpolation points = %d", inter_number))
    % fullFilePath = sprintf("./t4_%d.png", inter_number);
    % exportgraphics(gcf, fullFilePath, 'Resolution', 300);
end

function res = interplot_Newton(f_list, x_list)
    % f_list(i) stores the function value at the i-th interpolation point. 
    % x_list stand for the list of interpolation points.(i.e. x_1, x_2, ..., x_n)
    n = length(x_list);
    a_list = zeros(n, 1); a_list(1) = f_list(1);
    % a_list(i) stores f[x_1, ..., x_i]
    for i = 1:n - 1
    % each turn, we caculate the i+1-th order differential
    % and at this time, only the front n-i part stores valid data.
        temp = (f_list(2) - f_list(1)) / (x_list(1 + i) - x_list(1));
        a_list(i + 1) = temp;
        for j = 2:n - i
            f_list(j - 1) = temp;
            temp = (f_list(j + 1) - f_list(j)) / (x_list(j + i) - x_list(j));
        end
        f_list(n - i) = temp;
    end
    % the same as: $\sum_{i=0}^{n-1} (f[x_1, ..., x_{i+1}]\prod_{j=1}^{i} (t - x_j))$
    res = @(t) sum(arrayfun(@(i) a_list(i) * prod(t - x_list(1:i - 1)), 1:n));
end

% function res = interplot_Lagrange(f_list, x_list)
%     n = length(x_list);
%     res = @(t) sum(arrayfun(@(i) f_list(i) * prod((t - x_list([1:i - 1 i + 1:end])) ./ ...
%         (x_list(i) - x_list([1:i - 1 i + 1:end]))), 1:n));
% end
