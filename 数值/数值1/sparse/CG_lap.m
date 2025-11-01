clc;clear;
n = 100;
[F, T] = gener(n);

A = kron(eye(n - 1), T) + kron(T', eye(n - 1));
f = F(:);
u=CG(A, f, zeros((n -1)^2, 1), 1e-16, 1000);
U = vec2mat(u);
Laplace(U);
title('Laplace operator(CG method)');

function u=CG(A, f, u0, tol,maxiter)
    u=u0;k=1;r=f-A*u0;
    lu=dot(r,r);p=r;lu_old=lu;
    tol=tol*norm(f);err_list=zeros(1,maxiter);
    err_list(1)=sqrt(lu);
    while (err_list(k)>tol && k<=maxiter)
        if(k>1)
            p=r+lu/lu_old*p;
        end
        w=A*p;
        alpha=lu/dot(p,w);
        u=u+alpha*p;
        r=r-alpha*w;
        lu_old=lu;
        lu=dot(r,r);
        k=k+1;
        err_list(k)=sqrt(lu);
    end
    if(k==maxiter)
        disp('Maximum number of iterations reached');   
    else
        figure;
        plot(log10(err_list(1:k)));
        xlabel('Number of iterations');
        ylabel('log-sacle of Error');
        title("Conjugate Gradient Method's convergence history");
        fprintf('Number of iterations: %d\n',k);
    end
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
