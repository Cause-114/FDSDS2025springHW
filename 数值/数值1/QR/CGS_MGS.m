function [Q, R] = CGS(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        % for j = 1:i - 1
        %     R(j, i) = dot(A(1:end, j), A(1:end, i));
        % end
        R(1:i - 1, i) = A(:, 1:i - 1)' * A(:, i);
        A(:, i) = A(:, i) - A(:, 1:i - 1) * R(1:i - 1, i);
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
    end
    Q = A;
end

function [Q, R] = MGS(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    Q = A;
end

function [Q, R] = CGS2(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        % for j = 1:i - 1
        %     R(j, i) = dot(A(:, j), A(:, i));
        % end
        R(1:i - 1, i) = A(:, 1:i - 1)' * A(:, i);
        A(:, i) = A(:, i) - A(:, 1:i - 1) * R(1:i - 1, i);
        % for j = 1:i - 1
        %     A(1:end, i) = A(1:end, i) - dot(A(1:end, j), A(1:end, i)) * A(1:end, j);
        % end
        temp = A(:, 1:i - 1)' * A(:, i);
        R(1:i - 1, i) = R(1:i - 1, i) + temp;
        A(:, i) = A(:, i) - A(:, 1:i - 1) * temp;
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
    end
    Q = A;
end

function [Q, R] = MGS2(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        % R(i, i) = norm(A(1:end, i));
        % A(1:end, i) = A(1:end, i) / R(i, i);
        for j = 1:i - 1
            temp = dot(A(:, j), A(:, i));
            R(j, i) = R(j, i) + temp;
            A(:, i) = A(:, i) - A(:, j) * temp;
        end
        % A(1:end, i) = A(1:end, i) / norm(A(1:end, i));
        R(i, i) = norm(A(:, i));
        A(1:end, i) = A(1:end, i) / R(i, i);
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    Q = A;
end

% You can change the size of matrix here;
M = 40;
N = 20;
A1 = randi([1, M ^ 2], M, N);
A2 = randi([1, M ^ 2], M, N);
A = complex(A1, A2);

figure('Name', 'good vesion1')
[Q, R] = CGS(A);
subplot(2, 2, 1);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with CGS method")
colorbar

subplot(2, 2, 2);
imagesc(abs(Q * R - A))
title("abs of QR-A with CGS method")
colorbar

[Q, R] = MGS(A);
subplot(2, 2, 3);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with MGS method")
colorbar

subplot(2, 2, 4);
imagesc(abs(Q * R - A))
title("abs of QR-A with MGS method")
colorbar

figure('Name', 'good vesion2')
[Q, R] = CGS2(A);
subplot(2, 2, 1);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with CGS2 method")
colorbar

subplot(2, 2, 2);
imagesc(abs(Q * R - A))
title("abs of QR-A with CGS2 method")
colorbar

[Q, R] = MGS2(A);
subplot(2, 2, 3);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with MGS2 method")
colorbar

subplot(2, 2, 4);
imagesc(abs(Q * R - A))
title("abs of QR-A with MGS2 method")
colorbar

sigma = zeros(M, N);
sigma(1, 1) = 4 ^ - (N / 2);

for i = 2:N
    sigma(i, i) = sigma(i - 1, i - 1) * 4;
end

U = orth(randn(M) + 1i * randn(M));
V = orth(randn(N) + 1i * randn(N));
U = U * exp(1i * angle(det(U)));
V = V * exp(1i * angle(det(V)));
A = U * sigma * V';

figure('Name', 'bad vesion1')
[Q, R] = CGS(A);
subplot(2, 2, 1);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with CGS method")
colorbar

subplot(2, 2, 2);
imagesc(abs(Q * R - A))
title("abs of QR-A with CGS method")
colorbar

[Q, R] = MGS(A);
subplot(2, 2, 3);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with MGS method")
colorbar

subplot(2, 2, 4);
imagesc(abs(Q * R - A))
title("abs of QR-A with MGS method")
colorbar

figure('Name', 'bad vesion2')
[Q, R] = CGS2(A);
subplot(2, 2, 1);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with CGS2 method")
colorbar

subplot(2, 2, 2);
imagesc(abs(Q * R - A))
title("abs of QR-A with CGS2 method")
colorbar

[Q, R] = MGS2(A);
subplot(2, 2, 3);
imagesc(abs(Q' * Q - eye(N)))
title("abs of Q'*Q-I with MGS2 method")
colorbar

subplot(2, 2, 4);
imagesc(abs(Q * R - A))
title("abs of QR-A with MGS2 method")
colorbar
