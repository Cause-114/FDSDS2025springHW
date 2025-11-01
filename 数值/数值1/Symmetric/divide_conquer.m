n = 6;
z = rand(n, 1) .* 2 + 0.001;
dia = [1; 2; 3; 4; 5; 6];
dia = dia + sort(rand(n, 1));
A = diag(dia) + z * z';

% here we generate distinct and sorted diagonal elements for the matrix A.
% And they are quite far from each other, which makes the pucture more beautiful.
[~, V] = eig(A);
v = sort(diag(V));
% disp(v);
% disp(dia);

% In the following code, we plot the result part by part,in order to avoid the Inf value.
figure;
x = linspace(dia(1) - 1, dia(1) - 0.0001, 1000);
y = f(x, z, dia);
plot(x, y, 'b');
hold on;

for i = 2:n
    l = v(i - 1) - 0.999 * (v(i - 1) - dia(i - 1));
    r = v(i - 1) + 0.999 * (dia(i) - v(i - 1));
    x = [linspace(l, v(i - 1), 500), linspace(v(i - 1), r, 500)];
    y = f(x, z, dia);
    plot(x, y, 'b');
    hold on;
end

l = v(n) - 0.999 * (v(n) - dia(n));
r = v(n) + 1;
x = [linspace(l, v(n), 500), linspace(v(n), r, 500)];
y = f(x, z, dia);
plot(x, y, 'b', 'DisplayName', 'function');
hold on;

plot(v, zeros(n, 1), '.r', 'DisplayName', 'Real Eigenvalues');
legend ("fuction", "Real Eigenvalues");

function y = f(lambda, z, dia)
    y = zeros(length(lambda), 1);
    for k = 1:length(lambda)
        y(k) = 1 + z' * (z ./ (dia - lambda(k)));
    end
end
