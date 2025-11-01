function [Q, R] = cholesky(B)
    A = B' * B;
    n = size(A, 1);
    for i = 1:n
        A(i, i) = sqrt(A(i, i));
        A(i + 1:n, i) = A(i + 1:n, i) / conj(A(i, i));
        for j = i + 1:n
            A(j:n, j) = A(j:n, j) - A(j:n, i) * conj(A(j, i));
        end
    end
    R = zeros(n, n);
    for i = 1:n
        R(i, i:n) = A(i:n, i)';
    end
    Q = B / R;
end

function [Q, R] = Household(A)
    n = size(A, 1);
    Q = eye(n,n);
    m=size(A,2);
    R = zeros(n, m);
    for i = 1:m
        sub_x = dot(A(i + 1:n, i), A(i + 1:n, i));
        if (sub_x~=0)
            norm_x = sqrt(A(i, i) * conj(A(i, i)) + sub_x);
            v=A(i:n, i);v(1)=v(1)-norm_x;
            beta=2/(dot(v,v));
            A(i:n,i:m)=A(i:n,i:m)-beta*v*(v'*A(i:n,i:m));
            R(i,:)=A(i,:);
            Q(:,i:n)=Q(:,i:n)-(Q(:,i:n)*beta*v)*v';
        end 
    end
end
N=100;
A1 = randi([1, 10], N, N);
A2=randi([1,10],N,N);
A=complex(A1,A2);
% [Q,R]=cholesky(A);
% disp(Q);
% disp(R);
% [Q,R]=Household(A);
% disp(Q);
% disp(R);
[Q,R]=cholesky(A);
figure
subplot(2,2,1);
imagesc(abs(Q * Q' - eye(N)))
title("abs of Q*Q'-I with cholesky method(ill version)")
colorbar
subplot(2,2,2);
imagesc(abs(Q * R - A))
title("abs of Q*R-A with cholesky method(ill version)")
colorbar

[Q, R] = Household(A);
subplot(2,2,3);
imagesc(abs(Q * Q' - eye(N)))
title("abs of Q*Q'-I with Householder method(ill version)")
colorbar
subplot(2,2,4);
imagesc(abs(Q * R - A))
title("abs of Q*R-A with Householder method(ill version)")
colorbar

sigma = eye(N);
sigma(1,1)=8.8e-16;
for i=2:100
    sigma(i,i)=sigma(i-1,i-1)*2;
end
% 生成两个随机的复数酉矩阵
U = orth(randn(N) + 1i * randn(N));
V = orth(randn(N) + 1i * randn(N));

% 确保U和V确实是酉矩阵
U = U * exp(1i * angle(det(U)));
V = V * exp(1i * angle(det(V)));

% 构造复矩阵A
A = U * sigma * V';
[Q, R] = cholesky(A);
figure
subplot(2,2,1);
imagesc(abs(Q * Q' - eye(N)))
title("abs of Q*Q'-I with cholesky method(ill version)")
colorbar
subplot(2,2,2);
imagesc(abs(Q * R - A))
title("abs of Q*R-A with cholesky method(ill version)")
colorbar

[Q, R] = Household(A);
subplot(2,2,3);
imagesc(abs(Q * Q' - eye(N)))
title("abs of Q*Q'-I with Householder method(ill version)")
colorbar
subplot(2,2,4);
imagesc(abs(Q * R - A))
title("abs of Q*R-A with Householder method(ill version)")
colorbar

% function [Q, R] = MGS(A)
%     m=size(A,2);
%     R = zeros(m, m);
%     for i = 1:m
%         R(i, i) = norm(A(:, i));
%         A(:, i) = A(:, i) / R(i, i);
%         for j = i + 1:m
%             R(i, j) = dot(A(:, i), A(:, j));
%             A(:, j) = A(:, j) - R(i, j) * A(:, i);
%         end
%     end
%     Q = A;
% end

% A=randi([1,10],15,5);
% A(1:5,1:5)=zeros(5,5);
% disp(A);
% [Q,R]=Household(A);
% disp(Q);
% disp(R);
% [Q,R]=MGS(A);   
% disp(Q);
% disp(R);