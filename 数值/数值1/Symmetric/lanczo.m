n=25;
iter_list=2:n;
los_orth =zeros(length(iter_list),1);
val=randi([1,n],n,1);
A=diag(val);
U=orth(rand(n,n));
A=U'*A*U;
figure;
for i=1:length(iter_list)
    m=iter_list(i);
    [Q, T] = lanczos(A, m);
    [T, ~] = upper_iter(T);
    plot(diag(T),m*ones(m,1),".");
    hold on;
    los_orth(i) = log(norm(Q'*Q-eye(m),inf));
end
plot(val,zeros(n,1),".");
hold off;
title('eig convergence history of Tk(with 0 itertime as the true eig)');
ylabel('iteration number(k)');
xlabel('eigenvalue');
figure
plot(iter_list,los_orth);
title('Orthogonality loss of Lanczos algorithm');
xlabel('iteration number(k)');
ylabel("Orthogonality loss(log-scale)(norm(Q'*Q-eye(m),inf))");

% A = rand(n, n);
% A = A' + A;
% m=30;
% [Q, T] = lanczos(A, m);
% [gam, V] = upper_iter(T);
% figure
% imagesc(abs(Q'*Q-eye(m)));
% colorbar;
% figure
% imagesc(abs(T*V-V*gam));
% colorbar;
% function dispp(u,m)
%     u=diag(u);
%     [~,idx]=sort(abs(u),"descend");
%     u=u(idx);
%     disp(u(1:m)')
% end
function [T,V] = upper_iter(T)
    m = size(T, 1);
    V=eye(m);
    while (true)
        for i = 1:m - 1
            if (abs(T(i + 1, i)) < (abs(T(i, i)) + abs(T(i + 1, i + 1))) * 1e-15)
                T(i + 1, i) = 0;T(i, i + 1) = 0;
            end
        end
        for i = m - 1:-1:1
            if (T(i + 1, i)~=0)
                break;
            end
            m = m - 1;
        end
        if (m < 2)
            break;
        end
        l = m - 1;
        for i = l - 1:-1:1
            if (T(i + 1, i) == 0)
                break;
            end
            l = i;
        end
        [T(l:m, l:m),GG]= shift_qr_iter(T(l:m, l:m));
        V(:, l:m) = V(:, l:m) * GG;
    end
end
function [c, s] = givens(a, b)
    if(b==0)
        c=1;s=0;
    else
        if(abs(b)>abs(a))
            tau=a/b;
            s=1/sqrt(1+tau^2);
            c=tau*s;
        else
            tau=b/a;
            c=1/sqrt(1+tau^2);
            s=tau*c;
        end 
    end
end

function [T,GG] = shift_qr_iter(T)
    n=size(T,1);
    GG=eye(n);
    d=(T(n-1,n-1)-T(n,n))/2;
    miu=T(n,n)-T(n,n-1)^2/(d+sign(d))*(sqrt(d^2+T(n,n-1)^2));
    x=T(1,1)-miu;z=T(2,1);
    for k=1:n-1
        [c,s]=givens(x,z);
        G=[c s;-s c];
        T(k:k+1,:)=G*T(k:k+1,:);
        T(:,k:k+1)=T(:,k:k+1)*G';
        GG(:,k:k+1)=GG(:,k:k+1)*G';
        if(k<n-1)
            x=T(k+1,k);z=T(k+2,k);
        end
    end
end

function [Q, T] = lanczos(A, m)
    n = size(A, 1);
    Q = zeros(n, m);
    Q(:, 1) = randn(n, 1); Q(:, 1) = Q(:, 1) / norm(Q(:, 1));
    T = zeros(m, m);
    u = A * Q(:, 1);
    beta= 0;
    for i = 1:m - 1
        T(i, i) = dot(Q(:, i), u);
        u = u - T(i, i) * Q(:, i);
        beta = norm(u);

        if beta < 1e-10
            break;
        end

        Q(:, i + 1) = u / beta;
        T(i + 1, i) = beta; T(i, i + 1) = beta';
        u = A * Q(:, i + 1) - T(i, i + 1) * Q(:, i);
    end

    if (beta > 1e-10)
        T(m, m) = dot(Q(:, m), u);
    else
        warning("the algorithm break due to small beta");
    end

end
