n = 1000; precision = 1e-10;
rng(123);
N = 1000;
t1 = 0.5*(1+(-1).^(0:N-1))'/sqrt(N/2);
t2 = 0.5*(1+(-1).^(1:N))'/sqrt(N/2);
U = [t1, t2];
S = [1000, 0; 0, 999];
A = U*S*U' + 0.01*rand(N);
% A = rand(n, n);
eigvals = power_meth(A, precision);
eigvals1 = invers_iter_Rayleigh(A, 1000);
eigvals2 = invers_iter(A, 1004);
max_eig = max(abs(eig(A)));
fprintf('The precise maximum eigenvalue is: %.10f.\n', max_eig);
fprintf('compared to the answer get by power method: %.10f\n', eigvals);
fprintf('compared to the answer get by Rayleigh quotient iteration: %.10f\n', eigvals1);
fprintf('compared to the answer get by inverse iteration: %.10f\n', eigvals2);

function eigvals = power_meth(A, precision)
    max_iter = 1; n = size(A, 1);
    b = ones(n, 1); b_plus = A * b;
    lamda_list = zeros(10000, 1);
    [~, index] = max(abs(b_plus));
    lamda_list(1) = b_plus(index);
    while norm(b_plus - lamda_list(max_iter) * b) >precision && max_iter<10000
        b = b_plus / lamda_list(max_iter);
        b_plus = A * b;
        [~, index] = max(abs(b_plus));
        max_iter = max_iter + 1;
        lamda_list(max_iter) = b_plus(index);
    end

    if (max_iter==10000)
        warning('Power method does not converge.');
    end
    lamda_list = lamda_list(1:max_iter);
    eigvals = lamda_list(max_iter);
    fprintf('When set the precision to be 1e%d,', log10(precision));
    fprintf('iteration ends after %d iterations.\n', max_iter);
    plot(1:max_iter, log10(abs(lamda_list-eigvals)));
    xlabel('Iteration');
    ylabel('Eigenvalue loss(log version)');
    title('Power Method');
    grid on;
end
function eigvals = invers_iter_Rayleigh(A, start_eig)
    % tic;
    max_iter = 1; precision = 1e-10;
    eigvals = start_eig; n = size(A, 1);
    b = ones(n, 1) / sqrt(n);
    b_plus = (A - eye(n) * eigvals) \ b;
    delta = 1 / dot(b_plus, b);
    % inverse of Rayleigh quotient,excactly (b'*b)/(b'*(A-eig*I)^-1*b);
    eig_list = zeros(100, 1);
    eig_list(1) = eigvals;

    while (abs(delta) > precision && max_iter < 100)
        eigvals = eigvals + delta;
        b = b_plus / norm(b_plus);
        b_plus = (A-eye(n)*eigvals) \ b;
        delta = 1 / dot(b, b_plus);
        % inverse of Rayleigh quotient,excactly (b'*b)/(b'*(A-eig*I)^-1*b);
        max_iter = max_iter + 1;
        eig_list(max_iter) = eigvals;
    end
    if max_iter == 100
        warning('Maximum number of iterations reached');
    end
    % fprintf("time spent: %f seconds by Rayleigh quotient iteration\n", toc);

    eig_list = eig_list(1:max_iter);
    figure;
    plot(1:max_iter, log10(abs(eig_list-eigvals)));
    xlabel('Iteration');
    ylabel('Eigenvalue loss(log version)');
    title('Rayleigh Quotient Iteration');
    fprintf('get answer after %d iterations by Rayleigh quotient iteration\n', max_iter - 1);
end

function eigvals = invers_iter(A, start_eig)
    % tic;
    max_iter = 1; precision = 1e-10;
    eigvals = start_eig; n = size(A, 1);
    eig_list = zeros(100, 1);
    eig_list(1) = eigvals;
    A = A - start_eig * eye(n); b = ones(n, 1);
    b_plus = A \ b;
    [~, index] = max(abs(b_plus));
    delta = 1 / b_plus(index);

    while (abs(delta + start_eig - eigvals) > precision && max_iter < 100)
        eigvals = start_eig + delta;
        max_iter = max_iter + 1;
        eig_list(max_iter) = eigvals;
        b = b_plus * delta;
        b_plus = A \ b;
        [~, index] = max(abs(b_plus));
        delta = 1 / b_plus(index);
    end
    if max_iter == 100
        warning('Maximum number of iterations reached');
    end
    % fprintf("time spent: %f seconds by inverse iteration\n", toc);

    eig_list = eig_list(1:max_iter);
    figure;
    plot(1:max_iter, log10(abs(eig_list-eigvals)));
    xlabel('Iteration');
    ylabel('Eigenvalue loss(log version)');
    title('Inverse Iteration');
    fprintf('get answer after %d iterations by inverse iteration\n', max_iter - 1);
end
