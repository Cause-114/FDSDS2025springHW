function ans_list = Guas_part_pivo(A, b_list)
    n = size(A, 1);
    p_list = zeros(n, 1);
    for i = 1:n - 1
        [~, max_idx] = max(abs(A(i:n, i)));
        max_idx = max_idx + i - 1;
        [A(max_idx, 1:n), A(i, 1:n)] = deal(A(i, 1:n), A(max_idx, 1:n));
        p_list(i) = max_idx;
        [b_list(i), b_list(p_list(i))] = deal(b_list(p_list(i)), b_list(i));
        A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
        A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - A(i + 1:n, i) * A(i, i + 1:n);
    end
    for i = 1:n - 1
        b_list(i + 1:n) = b_list(i + 1:n) - A(i + 1:n, i) * b_list(i);
    end
    for i = n:-1:1
        b_list(i) = b_list(i) / A(i, i);
        b_list(1:i - 1) = b_list(1:i - 1) - b_list(i) * A(1:i - 1, i);
    end
    ans_list = b_list;
end

function [A, b, x] = generate(log_kapa, m, n, t1, t2)
    U = orth(randn(m));
    V = orth(randn(n));
    sigma = zeros(m, n);
    t = 10 ^ (log_kapa / (n - 1));
    sigma(1, 1) = t ^ (-n / 2);
    for i = 2:n
        sigma(i, i) = sigma(i - 1, i - 1) * t;
    end
    A = U * sigma * V';
    x1 = randn(n, 1);
    x2 = randn(m - n, 1);
    x1 = x1 / norm(x1) * t1;
    x2 = x2 / norm(x2) * t2;
    b = U * [x1; x2];
    % for i = 1:n
    %     sigma(i, i) = 1 / sigma(i, i);
    % end
    % x = V * sigma' * U' * b;
    x=A\b;
end

function [x] = Cholesky_ls(B, b)
    A = B' * B;
    n = size(A, 1);
    for i = 1:n
        A(1:i - 1, i) = zeros(i - 1, 1);
        A(i, i) = sqrt(A(i, i));
        A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
        for j = i + 1:n
            A(j:n, j) = A(j:n, j) - A(j:n, i) * A(j, i);
        end
    end
    x = A' \ ((A \ B') * b);
end

function [x] = Householder_ls(A, b)
    n = size(A, 2);
    for i = 1:n
        sub_x = dot(A(i + 1:end, i), A(i + 1:end, i));
        if (sub_x)
            norm_x = sqrt(A(i, i) * conj(A(i, i)) + sub_x);
            v = A(i:end, i);
            if A(i, i) >= 0
                v(1) = -sub_x / (A(i, i) + norm_x);
            else
                v(1) = A(i, i) - norm_x;
            end
            d = 2 / (v(1) ^ 2 + sub_x);
            b(i:end) = b(i:end) - d * v * dot(v, b(i:end));
            A(i:end, i:end) = A(i:end, i:end) - d * v * (v' * A(i:end, i:end));
        end
    end
    x = A(1:n, :) \ b(1:n);
end

function [x] = MGS_ls(A, b)
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
    x = R \ A' * b;
end

function [x] = augumented(A, b)
    m = size(A, 1);
    n = size(A, 2);
    A = [eye(m), A; A', zeros(n, n)];
    b = [b; zeros(n, 1)];
    x = Guas_part_pivo(A,b);
    x = x(m + 1:end);
end

M = 40; N = 20;
log_kapa = 0:0.1:15;
err = zeros(151, 4);
err1 = zeros(151, 4);

for i = 1:151
    [A, b, xx] = generate(log_kapa(i), M, N, 10000, 0.0001);
    x = Cholesky_ls(A, b);
    err(i, 1) = norm(x - xx);
    x = Householder_ls(A, b);
    err(i, 2) = norm(x - xx);
    x = MGS_ls(A, b);
    err(i, 3) = norm(x - xx);
    x = augumented(A, b);
    err(i, 4) = norm(x - xx);
end

for i = 1:151
    [A, b, xx] = generate(log_kapa(i), M, N, 0.0001, 10000);
    x = Cholesky_ls(A, b);
    err1(i, 1) = norm(x - xx);
    x = Householder_ls(A, b);
    err1(i, 2) = norm(x - xx);
    x = MGS_ls(A, b);
    err1(i, 3) = norm(x - xx);
    x = augumented(A, b);
    err1(i, 4) = norm(x - xx);
end

figure("Name", "b close to Range(A)")
plot(log_kapa, log10(err(:, 1)), 'DisplayName', 'cholesky')
hold on
plot(log_kapa, log10(err(:, 2)), 'DisplayName', 'householder')
hold on
plot(log_kapa, log10(err(:, 3)), 'DisplayName', 'MGS')
hold on
plot(log_kapa, log10(err(:, 4)), 'DisplayName', 'augumented')
legend show
xlabel("kapa of problem")
ylabel("norm of x-x*")
figure("Name", "b far from Range(A)")
plot(log_kapa, log10(err1(:, 1)), 'DisplayName', 'cholesky')
hold on
plot(log_kapa, log10(err1(:, 2)), 'DisplayName', 'householder')
hold on
plot(log_kapa, log10(err1(:, 3)), 'DisplayName', 'MGS')
hold on
plot(log_kapa, log10(err1(:, 4)), 'DisplayName', 'augumented')
legend show
xlabel("kapa of problem")
ylabel("norm of x-x*")
