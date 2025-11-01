Fro = zeros(99, 2);

for n = 2:100
    A = rand(n, n);
    A = A + A';
    [vec, val] = cyclic_jacb(A, 0);
    [~, v] = eig(A);
    Fro(n - 1, 1) = norm(A * vec - vec * diag(val), "fro");
    v = sort(diag(v));
    val = sort(val);
    Fro(n - 1, 2) = norm(v - val);
end

figure;
subplot(2, 1, 1);
plot(log10(Fro(:, 1)));
xlabel('n');
ylabel('Frobenius norm of A*vec-vec*val');
subplot(2, 1, 2);
plot(log10(Fro(:, 2)));
xlabel('n');
ylabel('error between eig and jacobi');

n=10;
A = rand(n, n);
A = A + A';
[~, ~] = classic_jacb(A, 1);

function f = E(A)
    f = norm(A, "fro") - norm(diag(A));
end

function [idx, idy] = find_max(A)
    n = size(A, 1);
    val = A(2, 1);
    idx = 2; idy = 1;
    for i = 1:n
        [m_val, m_id] = max(abs(A(i + 1:end, i)));
        if (m_val > val)
            val = m_val;
            idx = i + m_id;
            idy = i;
        end
    end
end

function [eig_vec, eig_val] = classic_jacb(A,plot_statu)
    n = size(A, 1);
    eig_vec = eye(n);
    [x, y] = find_max(A);
    cnt = 0;
    convergence = zeros(5000, 1);
    while (abs(A(x, y)) > 1e-14 && cnt < 5000)
        tao = (2 * A(x, y)) / (A(x, x) - A(y, y));
        c2 = 1 / sqrt(1 + tao ^ 2);
        c = sqrt((c2 + 1) / 2);
        s = -c2 * tao / (2 * c);
        % here cos(2a) > 0 and cos(a) > 0, so abs(a)<pi/4;
        % bellow is the wrong version,which make a > pi/4.

        % c2 = -1 / sqrt(1 + tao ^ 2);
        % c = sqrt((c2 + 1) / 2);
        % s = -c2 * tao / (2 * c);
        R = [c, s; -s c];
        A([x, y], :) = R' * A([x, y], :);
        A(:, [x, y]) = A(:, [x, y]) * R;
        eig_vec(:, [x, y]) = eig_vec(:, [x, y]) * R;
        [x, y] = find_max(A);        
        cnt = cnt + 1;
        convergence(cnt) = E(A);
    end

    if (cnt == 5000)
        warning('Jacobi iteration times exceed 1000, the result may not be accurate');
    end
    convergence = convergence(1:cnt);
    if plot_statu
        figure;
        plot(1:cnt,log10(convergence));
        xlabel('Iteration times');
        ylabel('Error');
        title('Convergence of Jacobi method');
    end
    fprintf('Jacobi iteration times: %d\n', cnt);
    eig_val = diag(A);
end

function [eig_vec, eig_val] = cyclic_jacb(A, plot_statu)
    n = size(A, 1);
    eig_vec = eye(n);
    cnt = 0;
    sigma = 1000;
    delta = E(A);
    convergence = zeros(100, 1);

    while (delta > 1e-15 && cnt < 100)
        flag = 1;
        for x = 1:n - 1
            for y = x + 1:n
                if (abs(A(x, y)) > delta)
                    flag = 0;
                    tao = (2 * A(x, y)) / (A(x, x) - A(y, y));
                    c2 = 1 / sqrt(1 + tao ^ 2);
                    c = sqrt((c2 + 1) / 2);
                    s = -c2 * tao / (2 * c);
                    R = [c, s; -s c];
                    A([x, y], :) = R' * A([x, y], :);
                    A(:, [x, y]) = A(:, [x, y]) * R;
                    eig_vec(:, [x, y]) = eig_vec(:, [x, y]) * R;
                end
            end
        end
        cnt = cnt + 1;
        convergence(cnt) = E(A);
        if (flag)
            delta = delta / sigma;
        end
    end
    if (cnt >= 100)
        warning('Jacobi iteration times exceed 100, the result may not be accurate');
    end
    convergence = convergence(1:cnt);
    if plot_statu
        figure;
        plot(1:cnt,log10(convergence));
        xlabel('Iteration times');
        ylabel('Error');
        title('Convergence of Jacobi method');
    end
    eig_val = diag(A);
end
