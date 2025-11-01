time_spend = zeros(2, 1000);
n_level = 1:1000;

for n = 1:1000
    A = randi([1, 10 * n], n, n);
    AA = A;
    p_list = zeros(1, n);
    q_list = zeros(1, n);
    % disp(A);
    % partial pivoting
    tic;
    [~,~,~]=lu(A);
    % for i = 1:n
    %     [max_val, max_idx] = max(abs(A(i:n, i)));
    %     max_idx = max_idx + i - 1;

    %     if A(max_idx, i) == 0
    %         disp("singular matrixs!")
    %         break;
    %     end

    %     [A(max_idx, 1:n), A(i, 1:n)] = deal(A(i, 1:n), A(max_idx, 1:n));
    %     p_list(i) = max_idx;
    %     A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
    %     A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - A(i + 1:n, i) * A(i, i + 1:n);
    % end

    % L = eye(n);
    % U = zeros(n, n);

    % for i = 1:n
    %     U(i, i:n) = A(i, i:n);
    %     L(i + 1:n, i) = A(i + 1:n, i);
    % end

    time_spend(1, n) = toc;

    A = AA;
    % complete pivoting
    % tic;
    % [L,U,P,Q]=lu(A);
    % for i = 1:n
    %     B = A(i:n, i:n);
    %     [max_val, idx] = max(abs(B(:)));
    %     [x, y] = ind2sub(size(B), idx);
    %     max_idx = [x + i - 1, y + i - 1];

    %     if A(max_idx(1), max_idx(2)) == 0
    %         disp("singular matrixs!")
    %         break;
    %     end

    %     [A(max_idx(1), 1:n), A(i, 1:n)] = deal(A(i, 1:n), A(max_idx(1), 1:n));
    %     [A(1:n, max_idx(2)), A(1:n, i)] = deal(A(1:n, i), A(1:n, max_idx(2)));
    %     q_list(i) = max_idx(2);
    %     p_list(i) = max_idx(1);
    %     A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
    %     A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - A(i + 1:n, i) * A(i, i + 1:n);
    % end

    % L = eye(n);
    % U = zeros(n, n);

    % for i = 1:n
    %     U(i, i:n) = A(i, i:n);
    %     L(i + 1:n, i) = A(i + 1:n, i);
    % end

    % time_spend(2, n) = toc;
    disp(n);

end

plot(log(n_level(7:1000)), log(time_spend(1, 7:1000)), 'r', 'DisplayName', 'particial pivoting');
hold on;
% plot(log(n_level(7:1000)), log(time_spend(2, 7:1000)), 'b', 'DisplayName', 'complete pivoting');
legend show;
title('log-log scale table of time and LU factorization(using well_developed standard function)')
xlabel('size of matrix')
ylabel('time spend')
% disp(L);
% disp(U);
% disp(L * U);
% disp(q_list);
% disp(p_list);
