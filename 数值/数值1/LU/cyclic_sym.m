function ans_list = solve_triangle(A, b_list)
    step = 1; n = size(A, 1);

    while step < n

        for i = 1:step * 2:n - 2 * step
            times1 = A(i, i + step) / A(i + step, i + step);
            A(i, [i, i + step, i + 2 * step]) = A(i, [i, i + step, i + 2 * step]) - times1 * A(i + step, [i, i + step, i + 2 * step]);
            b_list(i) = b_list(i) - times1 * b_list(i + step);
            times1 = A(i + 2 * step, i + step) / A(i + step, i + step);
            A(i + 2 * step, [i, i + step, i + 2 * step]) = A(i + 2 * step, [i, i + step, i + 2 * step]) - times1 * A(i + step, [i, i + step, i + 2 * step]);
            b_list(i + 2 * step) = b_list(i + 2 * step) - times1 * b_list(i + step);
        end

        res = mod(n - 1, 2 * step);

        if (res >= step)
            res = n - res;
            times1 = A(res, res + step) / A(res + step, res + step);
            b_list(res) = b_list(res) - times1 * b_list(res + step);
            A(res, [res, res + step]) = A(res, [res, res + step]) - times1 * A(res + step, [res, res + step]);
        end

        step = step * 2;
    end

    b_list(1) = b_list(1) / A(1, 1);

    while step >= 1

        for i = step + 1:2 * step:n - step
            b_list(i) = (b_list(i) - b_list(i - step) * A(i, i - step) - b_list(i + step) * A(i, i + step)) / A(i, i);
        end

        res = mod(n - 1, 2 * step);

        if (res >= step)
            res = n - res;
            b_list(res + step) = (b_list(res + step) - b_list(res) * A(res + step, res)) / A(res + step, res + step);
        end

        step = step / 2;
    end

    ans_list = b_list;
end

n = 100;
A = 8 * eye(n);
b_list = 15 * ones(n, 1);
b_list(1) = 9; b_list(n) = 14;

for i = 1:n - 1
    A(i, i + 1) = 1;
    A(i + 1, i) = 6;
end

ans_list = solve_triangle(A, b_list);
disp(ans_list);

A = 6 * eye(n);
b_list = 15 * ones(n, 1);
b_list(1) = 7; b_list(n) = 14;

for i = 1:n - 1
    A(i, i + 1) = 1;
    A(i + 1, i) = 8;
end

ans_list = solve_triangle(A, b_list);
disp(ans_list);
