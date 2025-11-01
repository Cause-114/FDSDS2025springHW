function [Q, R] = CGS(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        for j = 1:i - 1
            R(j, i) = dot(A(:, j), A(:, i));
        end
        A(:, i) = A(:, i) - A(:, 1:i - 1) * R(1:i - 1, i);
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
    end
    Q = A;
end

function [Q, R] = MGS(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    Q = A;
end

function [Q, R] = CGS2(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        for j = 1:i - 1
            R(j, i) = dot(A(:, j), A(:, i));
        end
        A(:, i) = A(:, i) - A(:, 1:i - 1) * R(1:i - 1, i);
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
        for j = 1:i - 1
            A(:, i) = A(:, i) - dot(A(:, j), A(:, i)) * A(:, j);
        end
        A(:, i) = A(:, i) / norm(A(:, i));
    end
    Q = A;
end

function [Q, R] = MGS2(A)
    n = size(A, 2);
    R = zeros(n, n);
    for i = 1:n
        R(i, i) = norm(A(:, i));
        A(:, i) = A(:, i) / R(i, i);
        for j = 1:i - 1
            A(:, i) = A(:, i) - dot(A(:, j), A(:, i)) * A(:, j);
        end
        A(:, i) = A(:, i) / norm(A(:, i));
        for j = i + 1:n
            R(i, j) = dot(A(:, i), A(:, j));
            A(:, j) = A(:, j) - R(i, j) * A(:, i);
        end
    end
    Q = A;
end
function [Q, R] = Cho_QR(B)
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

function [Q, R] = Hou_QR(A)
    n = size(A, 2);
    d_list = zeros(1, n);
    m=size(A,1);
    Q = zeros(m,m);
    R = zeros(n, n);
    for i = 1:n
        sub_x = dot(A(i + 1:end, i), A(i + 1:end, i));
        if (sub_x)
            ei = A(i, i) / abs(A(i, i));
            norm_x = sqrt(A(i, i) * conj(A(i, i)) + sub_x) * ei;
            A(i, i) = -sub_x / (A(i, i) + norm_x) * ei * ei;
            d_list(i) = 2 *A(i, i) * conj(A(i, i))/ (A(i, i) * conj(A(i, i)) + sub_x);
            A(i:end, i) = A(i:end, i) / A(i, i);
            A(i:end, i + 1:end) = A(i:end, i + 1:end) - d_list(i) * A(i:end, i) * (A(i:end, i)' * A(i:end, i + 1:end));
            A(i, i) = norm_x;
            A(i,i:end)=A(i,i:end)/ei;
            Q(i,i)=ei;
        elseif(A(i,i)~=0)
            ei =A(i,i)/abs(A(i,i));
            A(i,i:end)=A(i,i:end)/ei;
            Q(i,i)=ei;    
        end 
    end
    for i = n:-1:1
        R(i, i:n) = A(i, i:n);
        v = zeros(m, 1);
        v(i) = 1;
        v(i + 1:end) = A(i + 1:end, i);
        Q = Q - d_list(i) * (v * ((v')*Q));
    end
    Q=Q(:,1:n);
end

log_kapa=0:0.1:15;
lose_QQ=zeros(151,6);
lose_A=zeros(151,6);
M=100;N=50;
% You can change size of matrix here;  
sigma = zeros(M, N);
U = orth(randn(M) + 1i * randn(M));
V = orth(randn(N) + 1i * randn(N));
U = U * exp(1i * angle(det(U)));
V = V * exp(1i * angle(det(V)));
for i=1:151
    t=10^(log_kapa(i)/(N-1));
    sigma(1,1)=1;
    for j=2:N
        sigma(j,j)=sigma(j-1,j-1)*t;
    end
    A = U * sigma * V';
    [Q,R]=CGS(A);
    lose_QQ(i,1)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,1)=norm(Q*R-A,'fro');
    [Q,R]=MGS(A);
    lose_QQ(i,2)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,2)=norm(Q*R-A,'fro');
    [Q,R]=CGS2(A);
    lose_QQ(i,3)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,3)=norm(Q*R-A,'fro');
    [Q,R]=MGS2(A);
    lose_QQ(i,4)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,4)=norm(Q*R-A,'fro');
    [Q,R]=Cho_QR(A);
    lose_QQ(i,5)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,5)=norm(Q*R-A,'fro');
    [Q,R]=Hou_QR(A);
    lose_QQ(i,6)=norm(Q'*Q-eye(N),'fro');
    lose_A(i,6)=norm(Q*R-A,'fro');
end
figure
plot(log_kapa,log10(lose_QQ(:,1)),'DisplayName','CGS');
hold on;
plot(log_kapa,log10(lose_QQ(:,2)),'DisplayName','MGS');
hold on;
plot(log_kapa,log10(lose_QQ(:,3)),'DisplayName','CGS2');
hold on;
plot(log_kapa,log10(lose_QQ(:,4)),'DisplayName','MGS2');
hold on;
plot(log_kapa,log10(lose_QQ(:,5)),'DisplayName','Cholesky');
hold on;
plot(log_kapa,log10(lose_QQ(:,6)),'DisplayName','Householder');
title("Frobenius norm of Q'Q-I of different ways");
xlabel("log kapa");
ylabel("Frobenius norm(log version)");
legend show;
figure
plot(log_kapa,log10(lose_A(:,1)),'DisplayName','CGS');
hold on;
plot(log_kapa,log10(lose_A(:,2)),'DisplayName','MGS');
hold on;
plot(log_kapa,log10(lose_A(:,3)),'DisplayName','CGS2');
hold on;
plot(log_kapa,log10(lose_A(:,4)),'DisplayName','MGS2');
hold on;
plot(log_kapa,log10(lose_A(:,5)),'DisplayName','Cholesky');
hold on;
plot(log_kapa,log10(lose_A(:,6)),'DisplayName','Householder');
title("Frobenius norm of QR-A of different ways");
xlabel("log kapa");
ylabel("Frobenius norm(log version)");
legend show;
