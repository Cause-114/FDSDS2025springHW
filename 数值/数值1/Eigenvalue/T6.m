lambda=3;
A=[lambda,1;0,lambda];
B=[lambda,1;0,-lambda];
eigvals = power_meth(A, 1e-6);
eigvals1 = power_meth(B, 1e-6);
fprintf('The eigenvalues of A and B are %f and %f.\n', eigvals, eigvals1);

function eigvals = power_meth(A, precision)
    max_iter = 1;
    % b=[1;0];
    b = rand(size(A, 1), 1);b=b/norm(b);
    b_plus = A * b;
    lamda_list = zeros(100, 1);
    [~, index] = max(abs(b_plus));
    lamda_list(1) = b_plus(index);

    while norm(b_plus - lamda_list(max_iter) * b) > precision && max_iter<100
        b = b_plus / lamda_list(max_iter);
        b_plus = A * b;
        [~, index] = max(abs(b_plus));
        max_iter = max_iter + 1;
        lamda_list(max_iter) = b_plus(index);
    end

    if (max_iter==100)
        warning('Power method does not converge.');
    end
    lamda_list = lamda_list(1:max_iter);
    eigvals = lamda_list(max_iter);
    fprintf('When set the precision to be 1e%d,', log10(precision));
    fprintf('iteration ends after %d iterations.\n', max_iter);
    figure;
    plot(1:max_iter, lamda_list);
    xlabel('Iteration');
    ylabel('Eigenvalue');
    title('Power Method');
    grid on;
end