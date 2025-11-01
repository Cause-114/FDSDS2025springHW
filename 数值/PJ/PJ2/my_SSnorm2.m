% function [U, V, S,cnt] = my_svd(A, eps)
function [norm_] = my_SSnorm2(A)
    if size(A, 1) < size(A, 2)
        A = A';
    end
    [delta, gamma] = biodiag(A);
    delta = svd_iter(delta, gamma, eps);
    norm_=max(abs(delta));
end

function [delta] = svd_iter(delta, gamma, eps)
    n = length(delta);
    while (true)
        nor_ = sqrt(dot(delta, delta) + dot(gamma, gamma)) * eps;
        for i = 1:n
            if (i ~= n && abs(gamma(i)) < eps * (abs(delta(i)) + abs(delta(i + 1))))
                gamma(i) = 0;
            end
            if (abs(delta(i)) < nor_)
                delta(i) = 0;
            end
        end
        for i = n - 1:-1:1
            if (gamma(i) ~= 0)
                if (delta(i) == 0)
                    [c, s] = Givens(delta(i + 1), gamma(i));
                    delta(i + 1) = c * delta(i + 1) + s * gamma(i);
                    gamma(i) = 0;
                else
                    break;
                end
            end
            n = i;
        end
        if (n == 1)
            break;
        end
        l = n - 1;
        for i = l - 1:-1:1
            if (gamma(i) == 0)
                break;
            end
            if (delta(i) == 0)
                for j = i + 1:n
                    [c, s] = Givens(delta(j), gamma(i));
                    delta(j) = c * delta(j) + s * gamma(i);
                    if (j ~= n)
                        gamma(i) = -s * gamma(j);
                        gamma(j) = c * gamma(j);
                    else
                        gamma(i) = 0;
                    end
                end
                break;
            end
            l = i;
        end
        [delta(l:n), gamma(l:n - 1)] = pro_iterate(delta(l:n), gamma(l:n - 1));
    end
end

function [c, s] = Givens(a, b)
    if b == 0
        c = 1; s = 0;
        return
    end
    if abs(b) > abs(a)
        tao = a / b;
        s = 1 / sqrt(1 + tao ^ 2); c = s * tao;
    else
        tao = b / a;
        c = 1 / sqrt(1 + tao ^ 2); s = c * tao;
    end
end

function [delta, gamma] = pro_iterate(delta, gamma)
    n = length(delta);
    alpha = delta(n) ^ 2 + gamma(n - 1) ^ 2;
    if (n > 2)
        de = (delta(n - 1) ^ 2 + gamma(n - 2) ^ 2 - alpha) / 2;
    else
        de = (delta(n - 1) ^ 2 - alpha) / 2;
    end
    beta = delta(n - 1) * gamma(n - 1);
    miu = alpha - beta ^ 2 / (de + sign(de) * sqrt(de ^ 2 + beta ^ 2));
    % Wilkinson shift
    y = delta(1) ^ 2 - miu; z = delta(1) * gamma(1);
    for k = 1:n - 1
        [c, s] = Givens(y, z);
        if (k ~= 1)
            gamma(k - 1) = c * y + s * z;
        end
        y = c * delta(k) + s * gamma(k);
        gamma(k) = -s * delta(k) + c * gamma(k);
        z = s * delta(k + 1); delta(k + 1) = c * delta(k + 1);
        [c, s] = Givens(y, z);
        delta(k) = c * y + s * z;
        if (k ~= n - 1)
            y = gamma(k) * c + delta(k + 1) * s;
            delta(k + 1) = -s * gamma(k) + c * delta(k + 1);
            z = s * gamma(k + 1); gamma(k + 1) = c * gamma(k + 1);
        else
            tmp = gamma(k);
            gamma(k) = c * gamma(k) + s * delta(k + 1);
            delta(k + 1) = -s * tmp + c * delta(k + 1);
        end
    end
end

function [delta, gamma] = biodiag(A)
    % assume that m is greater than n;
    m = size(A, 1); n = size(A, 2);
    for k = 1:n
        sig = dot(A(k + 1:m, k), A(k + 1:m, k));
        if sig ~= 0
            nor_ = sqrt(sig + A(k, k) ^ 2);
            if (A(k, k) < 0)
                A(k, k) = A(k, k) - nor_;
            else
                A(k, k) = -sig / (nor_ + A(k, k));
            end
            beta = 2 / (sig + A(k, k) ^ 2);
            A(k:m, k + 1:n) = A(k:m, k + 1:n) - beta * A(k:m, k) * (A(k:m, k)' * A(k:m, k + 1:n));
            A(k, k) = nor_;
        end
        if (k < n - 1)
            sig = dot(A(k, k + 2:n), A(k, k + 2:n));
            if sig ~= 0
                nor_ = sqrt(sig + A(k, k + 1) ^ 2);
                if (A(k, k + 1) < 0)
                    A(k, k + 1) = A(k, k + 1) - nor_;
                else
                    A(k, k + 1) = -sig / (nor_ + A(k, k + 1));
                end
                beta = 2 / (sig + A(k, k + 1) ^ 2);
                A(k + 1:m, k + 1:n) = A(k + 1:m, k + 1:n) - (A(k + 1:m, k + 1:n) * A(k, k + 1:n)') * beta * A(k, k + 1:n);
                A(k, k + 1) = nor_;
            end
        end
    end
    delta = diag(A); gamma = diag(A, 1);
end
