function ans_list = Guas_elimi_direc(A, b_list)
    n = size(A, 1);
    % AA的高斯消去(直接干)
    for i = 1:n - 1
        A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
        A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - A(i + 1:n, i) * A(i, i + 1:n);
    end

    % 求解Ly=b
    for i = 1:n - 1
        b_list(i + 1:n) = b_list(i + 1:n) - A(i + 1:n, i) * b_list(i);
    end

    % 求解Ux=y
    for i = n:-1:1
        b_list(i) = b_list(i) / A(i, i);
        b_list(1:i - 1) = b_list(1:i - 1) - b_list(i) * A(1:i - 1, i);
    end

    ans_list = b_list;
end

function ans_list = Guas_part_pivo(A, b_list)
    n = size(A, 1);
    p_list = zeros(n, 1);% 置换矩阵P
    % A的高斯消去(列选主元)
    for i = 1:n - 1
        [~, max_idx] = max(abs(A(i:n, i)));
        max_idx = max_idx + i - 1;
        [A(max_idx, 1:n), A(i, 1:n)] = deal(A(i, 1:n), A(max_idx, 1:n));
        % 列选主元，并做交换。

        p_list(i) = max_idx;
        [b_list(i), b_list(p_list(i))] = deal(b_list(p_list(i)), b_list(i));
        % 列选主元中，对b做交换。

        A(i + 1:n, i) = A(i + 1:n, i) / A(i, i);
        A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - A(i + 1:n, i) * A(i, i + 1:n);
    end

    % 求解Ly=b
    for i = 1:n - 1
        b_list(i + 1:n) = b_list(i + 1:n) - A(i + 1:n, i) * b_list(i);
    end

    % 求解Ux=y

    for i = n:-1:1
        b_list(i) = b_list(i) / A(i, i);
        b_list(1:i - 1) = b_list(1:i - 1) - b_list(i) * A(1:i - 1, i);
    end

    ans_list = b_list;
end

err_list = zeros(1, 99);
err_list1 = zeros(1, 99);

for n = 2:100
    A = 6 * eye(n);
    b_list = 15 * ones(n, 1);
    b_list(1) = 7; b_list(n) = 14;

    for i = 1:n - 1
        A(i, i + 1) = 1;
        A(i + 1, i) = 8;
    end

    % 初始化构造三对角矩阵A与b;

    ans_list = Guas_part_pivo(A, b_list);
    err_list(n - 1) = norm(ans_list - ones(n, 1));
    ans1_list = Guas_elimi_direc(A, b_list);
    err_list1(n - 1) = norm(ans1_list - ones(n, 1));
    % 实际求得的解与理想情况的差的范数。
    % if n == 100
    %     disp(ans1_list); % 直接消元
    %     disp(ans_list); % 列选主元
    % end

end

plot(2:100, err_list, "r*", "DisplayName", "error of Guassian elimination with partial pivoting");
hold on;
plot(2:100, err_list1, "b-", "DisplayName", "error of Guassian elimination without pivoting");
legend show;
title("error of Guassian elimination")
xlabel("size of matrix")
ylabel("norm of the difference between pratical result and ideal ones")
% 图示化
