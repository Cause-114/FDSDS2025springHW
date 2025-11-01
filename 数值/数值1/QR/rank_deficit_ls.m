% It's noticable that we only find the min_norm_result in this code.
% If you want to get all the result, remember to call
% N=null(A,'r'), and x+N*c (c is free) is the answer. 
M = 40; N = 20; rA = N/2;
log_kapa = 0:0.1:15;
err = zeros(151, 1);
err1 = zeros(151, 1);

for i = 1:151
    [A, b, xx] = generate(log_kapa(i), M, N, rA, 10000, 0.0001);
    x = MGS2_ls2(A, b, rA);
    err(i) = norm(x - xx);
    [A, b, xx] = generate(log_kapa(i), M, N, rA, 0.0001, 10000);
    x = MGS2_ls2(A, b, rA);
    err1(i) = norm(x - xx);
end

plot(log_kapa, log10(err(:,1)), 'DisplayName', 'MGS(b close to Range(A),rank deficient)')
hold on
plot(log_kapa, log10(err1(:,1)), 'DisplayName', 'MGS(b far from Range(A),rank deficient)')
legend show
xlabel("kapa of problem")
ylabel("norm of x-x*")

function [A, b, x] = generate(log_kapa, m, n, rA, t1, t2)
    U = orth(randn(m));
    V = orth(randn(n));
    sigma = zeros(m, n);
    t = 10 ^ (log_kapa / (rA - 1));
    sigma(1, 1) = t ^ (-rA / 2);
    for i = 2:rA
        sigma(i, i) = sigma(i - 1, i - 1) * t;
    end
    A = U * sigma * V';
    x1 = randn(rA, 1);
    x2 = randn(m - rA, 1);
    x1 = x1 / norm(x1) * t1;
    x2 = x2 / norm(x2) * t2;
    b = U * [x1; x2];
    % I enhance the accuracy of x this time by simply call A\b.
    x = A \ b;
end
% Due to problem of time, I only listed the MGS2 method to solve
% rank deficient least squares problem.
function [x] = MGS2_ls2(A, b, rA)
    n = size(A, 2); k = 1;
    R = zeros(rA, n); P = eye(n);
    for i = 1:rA
        for j = 1:i - 1
            delta = dot(A(:, j), A(:, i));
            A(:, j) = A(:, j) - delta * A(:, i);
            R(j, i) = R(j, i) + delta;
        end
        % reorthogonalization. 
        R(i, i) = norm(A(:, i));
        while R(i, i) <= 1e-5 && k + rA <= n
            [R(:, i), R(:, k + rA)] = deal(R(:, i), R(:, k + rA));
            [A(:, i), A(:, k + rA)] = deal(A(:, i), A(:, k + rA));
            [P(:, i), P(:, k + rA)] = deal(P(:, i), P(:, k + rA));
            % swap.
            k = k + 1;
            for j = 1:i - 1
                delta = dot(A(:, j), A(:, i));
                A(:, j) = A(:, j) - delta * A(:, i);
                R(j, i) = R(j, i) + delta;
            end
            % reorthogonalization.
            R(i, i) = norm(A(:, i));
        end
        % We make sure the A(:,1:i) is full rank in the loop above.
        A(:, i) = A(:, i) / R(i, i);
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    x = R \ (A(:, 1:rA)' * b);
    x = P * x;
end
