A=randi([1,10], 3, 3);
eig_=QReig(A);
[~,v]=eig(A);
v=diag(v);
fprintf("Eigenvalues of A by system:\n");
disp(v);
fprintf("Eigenvalues of A by QR decomposition:\n");
disp(eig_);

function eig_list = QReig(A)
    max_iter = 1;tol=1e-6;
    eig_list = diag(A);
    [Q, R] = MGS2(A);
    A = R * Q;
    eig_list_plus=diag(A);
    while(norm(eig_list_plus - eig_list) > tol && max_iter <100)
        eig_list = eig_list_plus;
        [Q, R] = MGS2(A);
        A = R * Q;
        eig_list_plus=diag(A);
        max_iter = max_iter + 1;
    end
    if(max_iter == 100)
        warning("Maximum number of iterations reached!");
    end
    fprintf("With tolerance 1e%d,Number of iterations: %d\n",log10(tol), max_iter);
end


function [Q, R] = MGS2(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
        for j = 1:i - 1
            A(:, i) = A(:, i) - dot(A(:, j), A(:, i)) * A(:, j);
        end
        A(:, i) = A(:, i) / norm(A(:, i));
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    Q = A;
end