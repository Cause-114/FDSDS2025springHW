n = 7;
runge = @(x) 1 / (1 + 25 * x .^ 2);
sample = linspace(-1, 1, n);
y = arrayfun(runge, sample);
plot_points = linspace(-3, 3, 1000);

figure;
plot(plot_points, arrayfun(runge, plot_points), '--');
hold on

plot(sample, y, 'o');
hold on

plot(plot_points, guassinter(y, sample, plot_points));
legend('Runge function', 'Sample points', 'Interpolation')
title(sprintf('Guassian Interpolation, %d points', n))
grid on

function fvalue = guassinter(y, sample, plot_points)
    n = length(sample);
    guass = @(x, m) exp(- n * (x - m) .^ 2);
    fvalue = zeros(1, length(plot_points));
    A = zeros(n, n);
    for i = 1:n
        hd = @(x) guass(x, sample(i));
        A(i, :) = arrayfun(hd, sample);
    end
    c = y / A;
    for i = 1:n
        hd = @(x) guass(x, sample(i));
        fvalue = fvalue + c(i) * arrayfun(hd, plot_points);
    end
end

%% The following code is by changing \sigma of Guassian function.
% n = 10;
% runge = @(x) 1 / (1 + 25 * x .^ 2);
% sample = linspace(-1, 0, n);
% y = arrayfun(runge, sample);
% plot_points = linspace(-3, 3, 1000);
% figure;
% plot(plot_points, arrayfun(runge, plot_points), '--')
% hold on
% plot(sample, y, 'o')
% hold on
% plot(plot_points, guassinter(y, sample, plot_points));
% legend('Runge function', 'Sample points', 'Interpolation')
% title(sprintf('Guassian Interpolation, %d points', n));
% grid on

% function fvalue = guassinter(y, sample, plot_points)
%     guass = @(x, bt) exp(- bt * x^ 2);
%     n = length(sample); beta = linspace(n-1,n,n);
%     fvalue = zeros(1, length(plot_points));
%     A = zeros(n, n);
%     for i = 1:n
%         hd = @(x) guass(x, beta(i));
%         A(i, :) = arrayfun(hd, sample);
%     end
%     c = y / A;
%     for i = 1:n
%         hd = @(x) guass(x, beta(i));
%         fvalue = fvalue + c(i) * arrayfun(hd, plot_points);
%     end
% end

