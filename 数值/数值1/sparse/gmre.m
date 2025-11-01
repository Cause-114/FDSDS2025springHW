rng default
n=100;tol=1e-6;
A = sprandn(n,n,0.1)+8*speye(n);
b = rand(n,1);x0 = zeros(n, 1);

x = GMRES(A, b, n, tol,x0);
fprintf("with norm of Ax-b:\n")
disp(norm(A * x - b))

A=(A+A')/2;
x = GMRES_sys(A, b,n,tol,x0);
fprintf("with norm of Ax-b:\n")
disp(norm(A * x - b))

function x = GMRES_sys(A, b, maxiter, tol,x0)
    err_list = zeros(maxiter, 1);
    n = size(A, 1); 
    r = b - A * x0;
    Q = zeros(n, maxiter + 1);
    H = zeros(maxiter + 1, maxiter);
    % lancoz过程：A*Q=Q*H。
    rj = [norm(r); zeros(maxiter, 1)];
    % rj 是一个列向量，记录了残差。
    beta = rj(1);
    Q(:, 1) = r / rj(1);
    QH = zeros(maxiter + 1, maxiter);
    RH=zeros(maxiter,maxiter);
    % QH 与RH 是指对H的QR分解。方便进行最小二乘求解。
    for i = 1:maxiter
        Q(:, i + 1) = A * Q(:, i);
        l=max(1,i-1);
        for j = l:i
            H(j, i) = dot(Q(:, j), Q(:, i + 1));
            Q(:, i + 1) = Q(:, i + 1) - H(j, i) * Q(:, j);
        end
        H(i + 1, i) = norm(Q(:, i + 1));
        QH(l:i + 1, i) = H(l:i + 1, i);
        l=max(l-1,1);
        for j = l:i - 1
            RH(j, i) = dot(QH(1:j + 1, j), QH(1:j + 1, i));
            QH(1:j + 1, i) = QH(1:j + 1, i) - QH(1:j + 1, j) * RH(j, i);
            % 对称情况下H(l:i+1,i)只可能在QH的l:i-1列有投影。
        end
        RH(i, i) = norm(QH(1:i + 1, i));
        QH(1:i + 1, i) = QH(1:i + 1, i) / RH(i, i);
        rj(1:i + 1) = rj(1:i + 1) - dot(QH(1:i + 1, i), rj(1:i + 1)) * QH(1:i + 1, i);
        % norm(Ax-b)=norm(r-A*Q*y)=norm(Q'*r-H*y)=norm(QH'*norm(r)*e1-RH*y).
        % 所以，只需计算e1除去QH列向量空间外的范数即可求解原问题残差，即rj(1:i+1)。
        err_list(i)=norm(rj(1:i+1));
        if (err_list(i) < tol)
            fprintf('GMRES_sys converged in %d iterations.\n', i);
            y=RH(1:i,1:i) \ (QH(1,1:i)'*beta);
            % 最小二乘求解。
            x = x0 + Q(:, 1:i) * y(1:i);
            figure;
            plot(log10(err_list(1:i)));
            title('Convergence of GMRES\_sys');
            xlabel('Iteration');
            ylabel('Error');
            return;
        end
        Q(:, i + 1) = Q(:, i + 1) / H(i + 1, i);
    end
    warning('MATLAB:gmres:NoConvergence', 'The algorithm did not converge.');
    x = x0 + Q(:, 1:maxiter) *(RH(1:maxiter,1:maxiter) \ (QH(1,1:maxiter)'*beta));

end
function x = GMRES(A, b, maxiter, tol,x0)
    err_list = zeros(maxiter, 1);
    n = size(A, 1); 
    r = b - A * x0;
    Q = zeros(n, maxiter + 1);
    H = zeros(maxiter + 1, maxiter);
    % arnoldi过程：A*Q=Q*H。
    rj = [norm(r); zeros(maxiter, 1)];
    beta = rj(1);
    Q(:, 1) = r / rj(1);
    % rj 是一个列向量，记录了残差。
    QH = zeros(maxiter + 1, maxiter);
    RH=zeros(maxiter,maxiter);
    % QH 与RH 是指对H的QR分解。方便进行最小二乘求解。
    for i = 1:maxiter
        Q(:, i + 1) = A * Q(:, i);
        for j = 1:i
            H(j, i) = dot(Q(:, j), Q(:, i + 1));
            Q(:, i + 1) = Q(:, i + 1) - H(j, i) * Q(:, j);
        end
        H(i + 1, i) = norm(Q(:, i + 1));
        QH(1:i + 1, i) = H(1:i + 1, i);
        for j = 1:i - 1
            RH(j, i) = dot(QH(1:j + 1, j), QH(1:j + 1, i));
            QH(1:j + 1, i) = QH(1:j + 1, i) - QH(1:j + 1, j) * RH(j, i);
        end
        RH(i, i) = norm(QH(1:i + 1, i));
        QH(1:i + 1, i) = QH(1:i + 1, i) / RH(i, i);
        rj(1:i + 1) = rj(1:i + 1) - dot(QH(1:i + 1, i), rj(1:i + 1)) * QH(1:i + 1, i);
        err_list(i)=norm(rj(1:i+1));
        % norm(Ax-b)=norm(r-A*Q*y)=norm(Q'*r-H*y)=norm(QH'*norm(r)*e1-RH*y).
        % 所以，只需计算e1除去QH列向量空间外的范数即可求解原问题残差，即rj(1:i+1)。
        if (err_list(i) < tol)
            fprintf('GMRES converged in %d iterations.\n', i);
            y=RH(1:i,1:i) \ (QH(1,1:i)'*beta);
            x = x0 + Q(:, 1:i) * y(1:i);
            figure;
            plot(log10(err_list(1:i)));
            title('Convergence of GMRES');
            xlabel('Iteration');
            ylabel('Error');
            return;
        end
        Q(:, i + 1) = Q(:, i + 1) / H(i + 1, i);
    end
    warning('MATLAB:gmres:NoConvergence', 'The algorithm did not converge.');
    x = x0 + Q(:, 1:maxiter) *(RH(1:maxiter,1:maxiter) \ (QH(1,1:maxiter)'*beta));
end



% function x = GMRES1(A, b, maxiter, tol,x0)
%     n = length(b);
%     r = b - A*x0;
%     beta = norm(r);
%     Q = zeros(n, maxiter); % 正交基
%     Q(:,1) = r / beta;
%     H = zeros(maxiter+1, maxiter); % 上 Hessenberg 矩阵
%     e1 = [1; zeros(maxiter, 1)]; % 单位向量
%     for i = 1:maxiter
%         q = A * Q(:,i);
%         for j = 1:i
%             H(j,i) = q' * Q(:,j);
%             q = q - H(j,i) * Q(:,j);
%         end
%         H(i+1,i) = norm(q);
%         if H(i+1,i) == 0
%             break;
%         end
%         Q(:,i+1) = q / H(i+1,i);

%         % 解决最小二乘问题
%         e1beta = beta * e1(1:i+1);
%         y = H(1:i+1, 1:i) \ e1beta;

%         % 更新解
%         x = x0 + Q(:,1:i) * y;

%         % 检查收敛
%         residual = norm(b - A*x);
%         if residual < tol
%             fprintf('GMRES1 converged in %d iterations.\n', i);
%             break;
%         end
%     end
% end