max_size = 100;
err_list = zeros(1, max_size-1);

for n = 2:max_size
    A = gener(n);
    [V, ~] = eig(A);
    [V1, ~] = my_eig_upper(A);
    err_list(n - 1) = norm(V1 - V) / norm(V);
end

plot(1:max_size-1, err_list);
xlabel('n');
ylabel('Relative error');
title('Relative error of eigvec for upper triangular matrix');

function [V, D] = my_eig_upper(A)
    n = size(A, 1);
    D = zeros(n, 1);
    V = eye(n);
    for i = 1:n
        D(i) = A(i, i);
        V(1:i - 1, i) = -A(1:i - 1, i);
        for j = i - 1:-1:1
            V(j, i) = V(j, i) / (A(j, j) - D(i));
            V(1:j - 1, i) = V(1:j - 1, i) - V(j, i) * A(1:j - 1, j);
        end
        V(1:i, i) = V(1:i, i) / norm(V(1:i, i));
    end
end

function A = gener(n) 
    A = diag(randperm(3 * n, n));
    for i = 2:n
        A(1:i - 1, i) = randi([1, 3 * n], i - 1, 1);
    end
end
