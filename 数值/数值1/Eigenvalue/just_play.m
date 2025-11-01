A = [1 0 0; 2 cos(pi / 3) -sin(pi / 3); 3 sin(pi / 3) cos(pi / 3)];
process(A);
% with complex eigenvalues;
A = [1 2; 3 4];
process(A);
% small scale matrix;
A = [1 2 3; 4 5 6; 7 8 9];
process(A);
% singular matrix;
A = [1, 0, 0, 0; 4, 1, 0, 0; 7, 8, 2, 0; 10, 11, 12, 3];
process(A);
% with repeated eigenvalues;
A=rand(100,100);
process(A);
% large scale matrix;

% A=[0 0 0 0 1
%     1 0 0 0 0
%     0 1 0 0 0
%     0 0 1 0 0
%     0 0 0 1 0];
% [Q,T]=my_schur(A);
% disp(Q);
% disp(T);
% size_list=50:4:200;
% iter_times=zeros(length(size_list), 1);
% for i = 1:length(size_list)
%     A = rand(size_list(i), size_list(i));
%     [Q, T,cnt] = my_schur(A);
%     iter_times(i) = cnt;
% end
% figure;
% plot(size_list, iter_times);
% title("Number of iterations required for my\_schur to converge");
% xlabel("Size of matrix");
% ylabel("Number of iterations");

function process(A)
    n = size(A, 1);
    [Q, T] = my_schur(A);
    if (is_schur(T))
        fprintf("congratulations!\n");
    end
    if (n >= 5)
        figure;
        imagesc(abs(Q * T * Q' - A));
        title("abs(Q*T*Q' - A)");
        colorbar;
        figure;
        imagesc(abs(Q * Q' - eye(n)));
        title("abs(Q*Q' - I)");
        colorbar;
    else
        fprintf("my_schur result:\n");
        disp(T);
        fprintf("compared with matlab-schur result:\n");
        T = schur(A, "real");
        disp(T);
    end
end

function comp = is_schur(T)
    n = size(T, 1);
    sub_diag=zeros(n-2,1);
    L_part=zeros(n-1,1);
    for i = 2:n
        if(i > 2)
            L_part(i-2) = norm(T(i, 1:i - 2),1);
            if (L_part(i-2) > 1e-5)
                fprintf("Matrix is not Schur!");
                fprintf("with norm of T(%d, 1:%d) = %f\n", i, i - 2, L_part(i-2));
                comp = false;
                return;
            end
        end
        if (abs(T(i, i - 1)) > 1e-5)
            if (abs(T(i, i) - T(i - 1, i - 1)) > 1e-5)
                fprintf("Matrix is not Schur!");
                fprintf("with 2x2 block T(%d, %d) = %f while T(%d, %d) = %f\n", i, i, T(i, i), i - 1, i - 1, T(i - 1, i - 1));
                comp = false;
                return;
            end
        else
            sub_diag(i - 1) = T(i, i-1);
        end
    end
    if(n>10)
        fprintf("norm of L_part:%f\n",norm(L_part,1));
        fprintf("norm of sub_diag:%f\n",norm(sub_diag,1));
    end
    comp = true;
end
