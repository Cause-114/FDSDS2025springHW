function comp = is_schur(T)
    n = size(T, 1);

    for i = 2:n

        if (i > 2 && norm(T(i, i - 2)) > 1e-5)
            fprintf("Matrix is not Schur!");
            fprintf("with norm of T(%d, 1:%d) = %f\n", i, i - 2, norm(T(i, 1:i - 2)));
            comp = false;
            return;
        end

        if (abs(T(i, i - 1)) > 1e-5)

            if (abs(T(i, i) - T(i - 1, i - 1)) > 1e-5)
                fprintf("Matrix is not Schur!,T(%d, %d)= %f\n", i, i - 1, T(i, i - 1));
                fprintf("with 2x2 block T(%d, %d) = %f while T(%d, %d) = %f\n", i, i, T(i, i), i - 1, i - 1, T(i - 1, i - 1));
                comp = false;
                return;
            end

        end

    end

    comp = true;
end
