N = 10; 
A1 = rand(N, N); A2 = randi(N, N);
A = complex(A1, A2); A = A + A';
% make it hermitian.
b = rand(N, 1);
[Q, H] = lanczo(A, b);
figure;
imagesc(abs(A - Q * H * Q'));
title("error of lanczos algorithm(abs(A - Q * H * Q')");
colorbar;
figure;
imagesc(abs(Q *Q'-eye(N)));
title("orthogonality of Q");
colorbar;

function [Q, H] = lanczo(A, p)
    n = size(A, 1);
    Q = zeros(n, n);
    H = zeros(n, n);
    Q(:, 1) = p / norm(p);
    p = A * Q(:, 1);
    for i = 1:n
        H(i, i) = dot(Q(:, i), p);
        p = p - Q(:, i) * H(i, i);
        beta = norm(p);
        if (abs(beta) < 1e-10)
            break;
        end
        % for j = 1:i
        %     p = p - Q(:, j) * dot(Q(:, j), p);
        % end
        % I write the loop above to reorthogonalize the vector p(make thre result more precise).
        % Without this step, result canbe defusing when matrix is large.
        % You can also skip it when the scale of test matrix is quite small.
        H(i + 1, i) = beta; H(i, i + 1) = beta';
        p = p / beta;
        Q(:, i + 1) = p;
        p = A * Q(:, i + 1) - Q(:, i) * H(i, i + 1);
    end
end
