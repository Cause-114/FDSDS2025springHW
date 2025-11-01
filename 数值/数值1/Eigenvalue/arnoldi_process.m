% You can change the size of the problem here.
N = 1000;
dim = 30;
ek = zeros(1, dim);
ek(dim) = 1;

A1 = randi([1, N ^ 3], N, N);
A2 = randi([1, N ^ 3], N, N);
A = complex(A1, A2); b = randi([1,N^2],N, 1);

figure("Name",'lose of orthogonality');
[Q, ~, ~, ~] = arnodi_CGS(A, b, dim);
subplot(2, 2, 1);
imagesc(abs(Q' * Q - eye(dim)));
title('CGS');
colorbar;
[Q, ~, ~, ~] = arnodi_MGS(A, b, dim);
subplot(2, 2, 2);
imagesc(abs(Q' * Q - eye(dim)));
title('MGS');
colorbar;
[Q, ~, ~, ~] = arnodi_CGS2(A, b, dim);
subplot(2, 2, 3);
imagesc(abs(Q' * Q - eye(dim)));
title('CGS2');
colorbar;
[Q, ~, ~, ~] = arnodi_MGS2(A, b, dim);
subplot(2, 2, 4);
imagesc(abs(Q' * Q - eye(dim)));
title('MGS2');
colorbar;

% In the following, you can test the accuracy of the methods.
% figure('Name','lose of arnoldi process')
% [Q, H, p, beta] = arnodi_CGS(A, b, dim);
% subplot(2, 2, 1);
% imagesc(abs(A * Q - Q * H - p * beta * ek));
% title('CGS');
% colorbar;
% [Q, H, p, beta] = arnodi_MGS(A, b, dim);
% subplot(2, 2, 2);
% imagesc(abs(A * Q - Q * H - p * beta * ek));
% title('MGS');
% colorbar;
% [Q, H, p, beta] = arnodi_CGS2(A, b, dim);
% subplot(2, 2, 3);
% imagesc(abs(A * Q - Q * H - p * beta * ek));
% title('CGS2');
% colorbar;
% [Q, H, p, beta] = arnodi_MGS2(A, b, dim);
% subplot(2, 2, 4);
% imagesc(abs(A * Q - Q * H - p * beta * ek));
% title('MGS2');
% colorbar;

function [Q, H, p, beta] = arnodi_CGS(A, b, dim)
    n = size(A, 1);
    Q = zeros(n, dim);
    H = zeros(dim, dim);
    Q(:, 1) = b / norm(b);
    for i = 2:dim + 1
        p = A * Q(:, i - 1);
        for j = 1:i - 1
            H(j, i - 1) = dot(Q(:, j), p);
        end
        p = p - Q(:, 1:i - 1) * H(1:i - 1, i - 1);
        beta = norm(p);
        p = p / beta;
        if i <= dim
            H(i, i - 1) = beta;
            Q(:, i) = p;
        end
    end
end

function [Q, H, p, beta] = arnodi_MGS(A, b, dim)
    n = size(A, 1);
    Q = zeros(n, dim);
    H = zeros(dim, dim);
    Q(:, 1) = b / norm(b);
    for i = 2:dim + 1
        p = A * Q(:, i - 1);
        for j = 1:i - 1
            H(j, i - 1) = dot(Q(:, j), p);
            p = p - Q(:, j) * H(j, i - 1);
        end
        beta = norm(p);
        p = p / beta;
        if i <= dim
            H(i, i - 1) = beta;
            Q(:, i) = p;
        end
    end
end

function [Q, H, p, beta] = arnodi_CGS2(A, b, dim)
    n = size(A, 1);
    Q = zeros(n, dim);
    H = zeros(dim, dim);
    Q(:, 1) = b / norm(b);
    for i = 2:dim + 1
        p = A * Q(:, i - 1);
        for j = 1:i - 1
            H(j, i - 1) = dot(Q(:, j), p);
        end
        p = p - Q(:, 1:i - 1) * H(1:i - 1, i - 1);
        h_delta = zeros(i - 1, 1);
        for j = 1:i - 1
            h_delta(j) = dot(Q(:, j), p);
            H(j, i - 1) = h_delta(j) + H(j, i - 1);
        end
        p = p - Q(:, 1:i - 1) * h_delta;
        beta = norm(p);
        p = p / beta;
        if i <= dim
            H(i, i - 1) = beta;
            Q(:, i) = p;
        end
    end
end

function [Q, H, p, beta] = arnodi_MGS2(A, b, dim)
    n = size(A, 1);
    Q = zeros(n, dim);
    H = zeros(dim, dim);
    Q(:, 1) = b / norm(b);
    for i = 2:dim + 1
        p = A * Q(:, i - 1);
        for j = 1:i - 1
            H(j, i - 1) = dot(Q(:, j), p);
            p = p - Q(:, j) * H(j, i - 1);
        end
        for j = 1:i - 1
            delta = dot(Q(:, j), p);
            p = p - Q(:, j) * delta;
            H(j, i - 1) = H(j, i - 1) + delta;
        end
        beta = norm(p);
        p = p / beta;
        if i <= dim
            H(i, i - 1) = beta;
            Q(:, i) = p;
        end
    end
end