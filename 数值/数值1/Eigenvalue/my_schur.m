function [Q, T] = my_schur(A)
% function [Q, T,iter_time] = my_schur(A)
    [Q, H] = house_hessen(A);
    [Q, H] = upper_iter(Q, H);
    % [Q, H,iter_time] = upper_iter(Q, H);
    [Q, T] = adjust2by2(Q, H);
end

function [Q, A] = house_hessen(A)
    n = size(A, 1);
    Q = eye(n);
    for i = 1:n - 2
        sigma = dot(A(i + 2:n, i), A(i + 2:n, i));
        if sigma ~= 0
            nor_ = sqrt(sigma + A(i + 1, i) ^ 2);
            if A(i + 1, i) < 0
                A(i + 1, i) = A(i + 1, i) - nor_;
            else
                A(i + 1, i) = -sigma / (A(i + 1, i) + nor_);
            end
            beta = 2 / (sigma + A(i + 1, i) ^ 2);
            Q(:, i + 1:n) = Q(:, i + 1:n) - beta * (Q(:, i + 1:n) * A(i + 1:n, i)) * A(i + 1:n, i)';
            A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - beta * A(i + 1:n, i) * (A(i + 1:n, i)' * A(i + 1:n, i + 1:n));
            A(:, i + 1:n) = A(:, i + 1:n) - beta * (A(:, i + 1:n) * A(i + 1:n, i)) * A(i + 1:n, i)';
            A(i + 2:n, i) = 0;
            A(i + 1, i) = nor_;
        end
    end
end

function [Q, H] = upper_iter(Q, H)
% function [Q, H,cnt] = upper_iter(Q, H)
    m = size(H, 1);
    % cnt=0;
    while (true)
        for i = 1:m - 1
            if (abs(H(i + 1, i)) < (abs(H(i, i)) + abs(H(i + 1, i + 1))) * 1e-15)
                H(i + 1, i) = 0;
            end
        end
        for i = m - 1:-1:2
            if (H(i + 1, i) && H(i, i - 1))
                break;
            end
            m = m - 1;
        end
        if (m < 3)
            break;
        end
        l = m - 2;
        for i = l - 1:-1:1
            if (H(i + 1, i) == 0)
                break;
            end
            l = i;
        end
        [H(l:m, l:m), H(1:l - 1, l:m), H(l:m, m + 1:end), Q(:, l:m)] = double_shift_qr(H(l:m, l:m), H(1:l - 1, l:m), H(l:m, m + 1:end), Q(:, l:m));
        % cnt=cnt+1;
    end
    % fprintf('Number of iterations: %d ,while n=%d\n', cnt, size(Q,1));
end

function [H, H1, H2, Q] = double_shift_qr(H, H1, H2, Q)
    n = size(H, 1);
    s = H(n, n) + H(n - 1, n - 1);
    t = H(n - 1, n - 1) * H(n, n) - H(n - 1, n) * H(n, n - 1);
    x = H(1, 1) * H(1, 1) + H(1, 2) * H(2, 1) - s * H(1, 1) + t;
    y = H(2, 1) * (H(1, 1) + H(2, 2) - s);
    z = H(2, 1) * H(3, 2);
    for k = 0:n - 3
        sig = y ^ 2 + z ^ 2;
        if sig ~= 0
            v = [x; y; z];
            if x < 0
                v(1) = x - sqrt(sig + x ^ 2);
            else
                v(1) = -sig / (x + sqrt(sig + x ^ 2));
            end
            beta = 2 / (v(1) ^ 2 + sig);
            Q(:, k + 1:k + 3) = Q(:, k + 1:k + 3) - beta * (Q(:, k + 1:k + 3) * v) * v';
            H1(:, k + 1:k + 3) = H1(:, k + 1:k + 3) - beta * (H1(:, k + 1:k + 3) * v) * v';
            H2(k + 1:k + 3, :) = H2(k + 1:k + 3, :) - beta * v * (v' * H2(k + 1:k + 3, :));
            l = max([1, k]);
            H(k + 1:k + 3, l:n) = H(k + 1:k + 3, l:n) - beta * v * (v' * H(k + 1:k + 3, l:n));
            r = min([k + 4, n]);
            H(1:r, k + 1:k + 3) = H(1:r, k + 1:k + 3) - beta * (H(1:r, k + 1:k + 3) * v) * v';
        end
        x = H(k + 2, k + 1); y = H(k + 3, k + 1);
        if (k ~= n - 3)
            z = H(k + 4, k + 1);
        end
    end
    if y ~= 0
        v = [x; y];
        if x < 0
            v(1) = x - sqrt(y ^ 2 + x ^ 2);
        else
            v(1) = -y ^ 2 / (x + sqrt(y ^ 2 + x ^ 2));
        end
        beta = 2 / (v(1) ^ 2 + y ^ 2);
        Q(:, n - 1:n) = Q(:, n - 1:n) - beta * (Q(:, n - 1:n) * v) * v';
        H1(:, n - 1:n) = H1(:, n - 1:n) - beta * (H1(:, n - 1:n) * v) * v';
        H2(n - 1:n, :) = H2(n - 1:n, :) - beta * v * (v' * H2(n - 1:n, :));
        H(n - 1:n, n - 2:n) = H(n - 1:n, n - 2:n) - beta * v * (v' * H(n - 1:n, n - 2:n));
        H(1:n, n - 1:n) = H(1:n, n - 1:n) - beta * (H(1:n, n - 1:n) * v) * v';
    end
end

function [Q, H] = adjust2by2(Q, H)
    n = size(H, 1);
    for i = 1:n - 1
        if (H(i + 1, i))
            tr = trace(H(i:i + 1, i:i + 1)); de = det(H(i:i + 1, i:i + 1));
            delt = tr ^ 2 - 4 * de;
            if (delt >= 0)
                lambda = tr / 2;
                if (H(i,i) > H(i + 1, i + 1))
                    lambda = lambda + sqrt(delt) / 2;
                else
                    lambda = lambda - sqrt(delt) / 2;
                end
                q1 = [lambda-H(i + 1, i + 1); H(i + 1, i)];
                q1 = q1 / norm(q1);
                % for that H(i+1,i)!=0,so q1 is not zero vector,it's save to divide its norm.
                Q_ = [q1(1), -q1(2); q1(2), q1(1)];
                H(1:i + 1, i:i + 1) = H(1:i + 1, i:i + 1) * Q_;
                H(i:i + 1, i:end) = Q_' * H(i:i + 1, i:end);
                Q(:, i:i + 1) = Q(:, i:i + 1) * Q_;
            else
                alpha = tr / 2;
                beta1 = (H(i + 1, i) - H(i, i + 1)) / 2;
                if (tr < 0)
                    beta1 = beta1 + sqrt(beta1 ^ 2 * 4 + delt) / 2;
                else
                    beta1 = beta1 - sqrt(beta1 ^ 2 * 4 + delt) / 2;
                end
                q1 = [-H(i, i + 1)-beta1; H(i, i)-alpha];
                if (norm(q1) < 1e-10)
                    q1 = [alpha-H(i + 1, i + 1); H(i + 1, i)-beta1];
                    if (norm(q1) < 1e-10)
                    %so the matrix is already a right one,no need to adjust.
                        continue;
                    end
                end
                q1 = q1 / norm(q1);
                Q_ = [q1(1), -q1(2); q1(2), q1(1)];
                H(1:i + 1, i:i + 1) = H(1:i + 1, i:i + 1) * Q_;
                H(i:i + 1, i:end) = Q_' * H(i:i + 1, i:end);
                Q(:, i:i + 1) = Q(:, i:i + 1) * Q_;
            end
        end
    end
end
