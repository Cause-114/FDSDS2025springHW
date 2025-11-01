x = [-1.0000, -1.0000, 1.0000, 1.0000, -0.7313, 0.5275, -0.0091, 0.3031];
y = [-1.0000, 1.0000, -1.0000, 1.0000, 0.6949, -0.4899, -0.1010, 0.5774];
% z = [1.6389, 0.5403, -0.9900, 0.1086, 0.9573, 0.8270, 1.6936, 1.3670];
guass=@(x,y) exp(-x^2-y^2);
z=zeros(1,length(x));
for i=1:length(x)
    z(i)=guass(x(i),y(i));
end
pt = complex(x, y);
xi = linspace(-1, 1, 100);
yi = linspace(-1, 1, 100);
[X, Y] = meshgrid(xi, yi);

% this one is working by using Newton interpolation with complex number.
figure;
Z = complex_2d(z, pt, xi, yi);
surf(X, Y, Z); xlabel('x');
ylabel('y'); title('Newton Interpolation');

% this one is working by using Lagrange interpolation 
% with reciprocal of distance as weight.
figure;
Z = lag_2d(z, pt, xi, yi);
surf(X, Y, Z); xlabel('x');
ylabel('y'); title('Lagrange Interpolation');

% this one is working by using delaunay triangulation.
figure;
tri = delaunay(x, y);
trisurf(tri, x, y, z, 'FaceColor', 'interp', 'EdgeColor', 'none');
xlabel('x'); ylabel('y');
title('Delaunay Triangulation Interpolation');

% this one is working by using cubic spline interpolation.
figure;
Z = griddata(x, y, z, X, Y, 'cubic');
surf(X, Y, Z); xlabel('x'); ylabel('y');
title('Cubic Spline Interpolation');

function res = complex_2d(f, x, xi, yi)
    n = length(x);
    for i = 1:n - 1
    % each loop rise one degree of differential.
    % and the first k-degree differential is f(k+1).
        tmp2 = f(i);
        % the first i f stores the previous first k-degree differential.
        % k=0,1,...i-1
        for j = 1:n - i
            tmp1 = (f(j + i) - f(i + j - 1)) / (x(j + i) - x(j));
            f(j + i - 1) = tmp2;
            tmp2 = tmp1;
        end
        f(n) = tmp2;
    end
    res = zeros(length(yi), length(xi));
    for i = 1:size(res, 1)
        for j = 1:size(res, 2)
            % calculate the interpolation value at complex number p.
            p = complex(xi(j), yi(i));
            prod_temp = 1;
            for k = 1:n - 1
                res(i, j) = res(i, j) + f(k) * prod_temp;
                prod_temp = prod_temp * (p - x(k));
            end
            res(i, j) = res(i, j) + f(n) * prod_temp;
        end
    end
    res = real(res);
end

function res = lag_2d(f, x, xi, yi)
    % You can change the rule of distance here.
    % I use the 2-norm here.
    mynorm = @ (x, p) (abs(real(x)) ^ p + abs(imag(x)) ^ p) ^ (1 / p);
    n = length(x);
    res = zeros(length(yi), length(xi));
    for i = 1:size(res, 1)
        for j = 1:size(res, 2)
            p = complex(xi(j), yi(i));
            % s stands for sum of reciprocal of distance.
            % flag aims to check if p is one of the data points.
            s = 0; flag = 0;
            for k = 1:n
                dis = mynorm(p - x(k), 3);
                if (dis == 0)
                % then it will become a disaster 
                % if we caculate reciprocal of distance.
                    res(i, j) = f(k);
                    flag = 1;
                    break;
                end
                res(i, j) = res(i, j) + f(k) / dis;
                s = s + 1 / dis;
            end
            if (flag == 0)
            % if no break, remember to regularize the result
            % by dividing it by sum of reciprocal of distance.
                res(i, j) = res(i, j) / s;
            end
        end
    end
end
