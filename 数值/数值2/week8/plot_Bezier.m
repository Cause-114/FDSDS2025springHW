% Write a program to plot the B´ezier curve with arbitrary set of control points.
% As it's better to visualize the curve in a 2D platform, we assume the input is 2D.

clc; clear; close all;
P=[0 0; 1 1; 2 0];
PLOT_BEZIER_CURVE(P); 
figure;
P = [0 0; 1 0; 1 1; 0 1]; % another example
PLOT_BEZIER_CURVE(P); 
function PLOT_BEZIER_CURVE(P)
    % P is a n by 2 matrix of control points, each row of
    % P is a control point (x,y)
    t = 0:0.01:1;
    x = zeros(size(t));
    y = zeros(size(t));
    for i = 1:length(t)
        [x(i), y(i)] = bezier(t(i), P);
    end
    plot(x, y, P(:, 1), P(:, 2), 'o');
    legend('Bezier Curve', 'Control Points');
end

function [x, y] = bezier(t, P)
    % P is a n by 2 matrix of control points, each row of
    % P is a control point (x,y), t is the parameter value.
    n = size(P, 1);
    for i = 1:n - 1
        for j = 1:n - i
            P(j, :) = (1 - t) * P(j, :) + t * P(j + 1, :);
        end
    end
    x = P(1, 1); y = P(1, 2);
end
