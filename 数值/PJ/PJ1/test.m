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

% size_list=50:200;
% zero_elements = zeros(length(size_list),1);
% for i = 1:length(size_list)
%     A = rand(size_list(i), size_list(i));
%     [Q,T]=my_schur(A);
%     [comp, zero_elements(i)] = is_schur(T);
% end
% figure;
% plot(size_list, log10(zero_elements));
% title("the sum of abs of zero elements in the left-bellow diagonal of T");
% xlabel("Size of matrix");
% ylabel("log10(abs of zero elements sum)");

% size_list=50:200;
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

% size_list=200:500;
% times=zeros(length(size_list), 1);
% for i = 1:length(size_list)
%     A = rand(size_list(i), size_list(i));
%     tic;
%     [Q, T] = my_schur(A);
%     times(i) = toc;
% end
% figure;
% plot(log10(size_list), log10(times));
% title("log10-sacle time required for my\_schur to work");
% xlabel("log10-Size of matrix");
% ylabel("log10-time(seconds)");

function process(A)
    n = size(A, 1);
    tic;
    [Q, T] = my_schur(A);
    fprintf("time spent: %f seconds\n", toc);
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
% function [comp, all_zero_sum] = is_schur(T)
    n = size(T, 1);all_zero_sum=0;epsilon=1e-10;
    for i = 2:n
        if (i > 2 && norm(T(i,1:i - 2),1) > epsilon)
            fprintf("Matrix is not Schur!");
            fprintf("with norm of T(%d, 1:%d) = %f\n", i, i - 2, norm(T(i, 1:i - 2)));
            comp = false;
            return;
        end
        all_zero_sum = all_zero_sum + norm(T(i, 1:i - 2), 1);
        if (abs(T(i, i - 1)) > epsilon)
            if (abs(T(i, i) - T(i - 1, i - 1)) > 1e-6)
                fprintf("Matrix is not Schur!");
                fprintf("with 2x2 block T(%d, %d) = %f while T(%d, %d) = %f\n", i, i, T(i, i), i - 1, i - 1, T(i - 1, i - 1));
                comp = false;
                return;
            end
        else
            all_zero_sum = all_zero_sum + abs(T(i, i - 1));
        end
    end
    if(n>5)
        fprintf("The elements of T are expected to be zero all sum up to 10^%f\n", log10(all_zero_sum));
        fprintf("with the size of T being %d\n", n);
    end
    comp = true;
end
