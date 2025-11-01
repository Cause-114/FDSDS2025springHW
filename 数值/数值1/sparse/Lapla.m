clc;clear;
n = 10;
[F, T] = gener(n);

A = kron(eye(n - 1), T) + kron(T', eye(n - 1));
f = F(:);
u=Guass_Seidel(A,f,zeros((n-1)^2,1),1e-12);
U = vec2mat(u);
Laplace(U);
title('Laplace operator(Guass-Seidel method)');

u=Jacobi(A,f,zeros((n-1)^2,1),1e-12);
U = vec2mat(u);
Laplace(U);
title('Laplace operator(Jacobi method)');



function x=Guass_Seidel(A,b,x0,tol)
    n=size(A,1);
    x=x0;
    err_list = zeros(1,1000);
    iter=1;
    for i=1:1000
        for j=1:n
            x(j)=(b(j)-A(j,1:j-1)*x(1:j-1)-A(j,j+1:n)*x(j+1:n))/A(j,j);
        end
        err_list(iter)=norm(A*x-b);
        if err_list(iter)<tol
            break
        end
        iter=iter+1;
    end
    if(iter==1000)
        warning('Guass-Seidel method does not converge');
    end
    figure;
    plot(log10(err_list(1:iter-1)));
    title('Convergence of Guass-Seidel method');
    xlabel('Number of iterations');
    ylabel('log10(error)');
end

function x=Jacobi(A,b,x0,tol)
    n=size(A,1);
    x=x0;
    err_list = zeros(1,1000);
    iter=1;
    for i=1:1000
        for j=1:n
            x(j)=(b(j)-A(j,1:j-1)*x0(1:j-1)-A(j,j+1:n)*x0(j+1:n))/A(j,j);
        end
        err_list(iter)=norm(A*x-b);
        if(err_list(iter)<tol)
            break
        end
        x0=x;
        iter=iter+1;
    end
    if(i==1000)
        warning('Jacobi method does not converge');
    end
    figure;
    plot(log10(err_list(1:iter-1)));
    title('Convergence of Jacobi method');
    xlabel('Number of iterations');
    ylabel('log10(error)');
end



function [F, T] = gener(n)
    F = zeros(n - 1, n - 1);
    F(1, :) = sin(pi * (1:n - 1) / n);
    T = 2 * eye(n - 1);
    for i = 1:n - 2
        T(i + 1, i) = -1;
        T(i, i + 1) = -1;
    end
end

function Laplace(U)
    n=size(U,1)-1;
    L = zeros(n - 1, n - 1);
    for i = 1:n - 1
        for j = 1:n - 1
            L(i, j) = (4 * U(i + 1, j + 1) - U(i, j + 1) - U(i + 2, j + 1) - U(i + 1, j) - U(i + 1, j + 2))*n^2;
        end
    end
    figure;
    imagesc(L);
    colorbar;
end

function U = vec2mat(u)
    n = sqrt(length(u)) + 1;
    U = zeros(n + 1, n + 1);
    U(1, 2:n) = sin(pi * (1:n - 1) / n);
    U(2:n, 2:n) = reshape(u, n - 1, n - 1);
end
