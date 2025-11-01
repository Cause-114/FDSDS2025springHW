n = 200;
A = hilb(n);
eigvals = invers_iter(A, 1);
eigvals2 = invers_iter_Rayleigh(A, 1);
[~, D] = eig(A);
eig_val = diag(D);
[~, index] = min(abs(eig_val - 1));
fprintf('The precise eigenvalue closest to 1 is: %.10f.\n', eig_val(index));
fprintf('Comparing to the answer by inverse iteration: %.10f\n', eigvals);
fprintf('Comparing to the answer by Rayleigh quotient iteration: %.10f\n', eigvals2);

function eigvals = invers_iter(A, start_eig)
    tic;
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
    fprintf("time spent: %f seconds by inverse iteration\n", toc);

    eig_list = eig_list(1:max_iter);
    figure;
    plot(1:max_iter, log10(abs(eig_list-eigvals)));
    xlabel('Iteration');
    ylabel('Eigenvalue loss(log version)');
    title('Inverse Iteration');
    fprintf('get answer after %d iterations by inverse iteration\n', max_iter - 1);
end

function eigvals = invers_iter_Rayleigh(A, start_eig)
    tic;
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
    fprintf("time spent: %f seconds by Rayleigh quotient iteration\n", toc);

    eig_list = eig_list(1:max_iter);
    figure;
    plot(1:max_iter, log10(abs(eig_list-eigvals)));
    xlabel('Iteration');
    ylabel('Eigenvalue loss(log version)');
    title('Rayleigh Quotient Iteration');
    fprintf('get answer after %d iterations by Rayleigh quotient iteration\n', max_iter - 1);
end