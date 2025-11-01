err1 = zeros(99, 2);
err2 = zeros(99, 2);

for n = 2:100
    A = complex(rand(n,n),rand(n, n));
    [Q, H] = Arnodi_Hessen(A);
    err1(n - 1, 1) = log10(norm(Q' * A * Q - H, 'fro'));
    err2(n - 1, 1) = log10(norm(eye(n) - Q * Q', 'fro'));
    [Q, H] = House_Hessen(A);
    err1(n - 1, 2) = log10(norm(Q' * A * Q - H , 'fro'));
    err2(n - 1, 2) = log10(norm(eye(n) - Q * Q', 'fro'));
end

figure;
subplot(2, 1, 1);
plot(err1(:, 1), 'DisplayName', 'Arnodi')
hold on;
plot(err1(:, 2), 'DisplayName', 'Householder')
legend show;
title("Frobenius norm of the Q'*A*Q-H");
xlabel('n');
ylabel('log10-error');
subplot(2, 1, 2);
plot(err2(:, 1), 'DisplayName', 'Arnodi')
hold on;
plot(err2(:, 2), 'DisplayName', 'Householder')
legend show;
title("Frobenius norm of the Q*Q'-I");
xlabel('n');
ylabel('log10-error');

function [Q, H] = Arnodi_Hessen(A)
    n = size(A, 1);
    Q = zeros(n, n);
    H = zeros(n, n);
    Q(:, 1) = rand(n, 1);
    Q(:, 1) = Q(:, 1) / norm(Q(:, 1));
    for i = 2:n
        Q(:, i) = A * Q(:, i - 1);
        for j = 1:i - 1
            H(j, i - 1) = dot(Q(:, j), Q(:, i));
            Q(:, i) = Q(:, i) - Q(:, j) * H(j, i - 1);
        end
        H(i, i -1) = norm(Q(:, i));
        Q(:, i) = Q(:, i) / H(i, i - 1);
    end
    p = A * Q(:, n);
    for i = 1:n
        H(i, n) = dot(Q(:, i), p);
        p = p - Q(:, i) * H(i, n);
    end
end

function [Q, H] = House_Hessen(A)
    n = size(A, 1);
    Q = eye(n);
    for i = 1:n - 2
        sigma = dot(A(i + 2:n, i), A(i + 2:n, i));
        if sigma ~= 0
            ei = A(i + 1, i) / abs(A(i + 1, i));
            nor_ = sqrt(sigma + A(i + 1, i) * A(i + 1, i)') * ei;
            A(i + 1, i) = -sigma / (A(i + 1, i) + nor_) * ei * ei;
            beta = 2 / (sigma + A(i + 1, i) * A(i + 1, i)');
            Q(:, i + 1:n) = Q(:, i + 1:n) - beta * (Q(:, i + 1:n) * A(i + 1:n, i)) * A(i + 1:n, i)';
            A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - beta * A(i + 1:n, i) * (A(i + 1:n, i)' * A(i + 1:n, i + 1:n));
            A(:, i + 1:n) = A(:, i + 1:n) - beta * (A(:, i + 1:n) * A(i + 1:n, i)) * A(i + 1:n, i)';
            A(i + 2:n, i) = 0;
            A(i + 1, i) = nor_;
        end
    end
    H = A;
end
