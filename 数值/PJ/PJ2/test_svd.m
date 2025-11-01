err_info = ["The sig is not sorted in descending order.\n", "The sig has negative values.\n", ...
                "U or V is not orthogonal.\n", "USV^T-A is not a matrix with zero mean.\n"];
pass_info = ["First 9 test cases(small size)!\n", "10~18 test cases(median size)!\n", ...
                "Finall test cases passed!(with different kapa)\n"];
testord = 1; stage = 1; flag = 0;
while stage <= 3
    while (testord <= 9 * stage)
        name = ['test', num2str(testord), '.mat'];
        load(['testcases\' name], 'A');
        [errord] = test(A, 1e-12);
        if errord > 0
            fprintf("Test case %d failed.\n", testord);
            fprintf(err_info(errord));
            flag = 1;
            break;
        end
        testord = testord + 1;
    end
    if (flag ~= 1)
        fprintf(pass_info(stage));
    else
        break;
    end
    stage = stage + 1;
end

function errord = test(A, tol)
    if nargin < 2
        tol = 1e-12;
    end
    errord = 0;
    [U, V, S] = my_svd(A, tol);
    m = size(A, 1); n = size(A, 2);
    for i = 2:min(m, n)
        if (S(i, i) > S(i - 1, i - 1))
            errord = 1;
            return;
        end
        if (S(i, i) < 0)
            errord = 2;
            return;
        end
    end
    tol2 = 1e3 * tol;
    if (norm(U' * U - eye(m)) > tol2 || norm(V * V' - eye(n)) > tol2)
        errord = 3;
        return;
    end
    if (norm(U * S * V' - A) > tol2)
        disp(norm(U * S * V' - A))
        errord = 4;
        return;
    end
end

%% bellow is the generated code for test cases
% ord = 1;
% for i = 1:3
%     for j = 1:3
%         A = rand(i, j);
%         name = ['test', num2str(ord), '.mat'];
%         save(['testcases\' name], 'A');
%         ord = ord + 1;
%     end
% end
% for i = 50:25:100
%     for j = 50:25:100
%         A = rand(i, j);
%         name = ['test', num2str(ord), '.mat'];
%         save(['testcases\' name], 'A');
%         ord = ord + 1;
%     end
% end
% for log_kapa = 1:9
%     m = randi([5, 50]); n = randi([5, 50]);
%     U = orth(randn(m)); V = orth(randn(n));
%     sigma = zeros(m, n); n = min(m, n);
%     t = 5 ^ (log_kapa / (n - 1));
%     sigma(1, 1) = t ^ (-n / 2);
%     for j = 2:n
%         sigma(j, j) = sigma(j - 1, j - 1) * t;
%     end
%     A = U * sigma * V';
%     name = ['test', num2str(ord), '.mat'];
%     save(['testcases\' name], 'A');
%     ord = ord + 1;
% end
%% test the function my_svd
% m=114;n=514;
% A=rand(m,n);
% [U, V, S] = my_svd(A, 1e-12);
% imagesc(U * S * V' - A);
% title('residual of SVD(USV^T-A)');
% colorbar;
% figure;
% subplot(1, 2, 1);
% imagesc(U * U' - eye(m));
% title('U*U^T-I');
% colorbar;
% subplot(1, 2, 2);
% imagesc(V * V' - eye(n));
% title('V*V^T-I');
% colorbar;

% size_n=50:250;
% iter_time=zeros(1,length(size_n));
% for i=1:length(size_n)
%     A=rand(size_n(i),size_n(i));
%     tic;
%     [~,~,~]=my_svd(A);
%     iter_time(i)=toc;
% end
% figure;
% plot(log10(size_n),log10(iter_time));
% xlabel('log10(size of n)');
% ylabel('log10(time)');
% title('Time required by my\_svd for different matrix n*n');

% size_n=2:100;
% iter_time=zeros(1,length(size_n));
% for i=1:length(size_n)
%     A=rand(size_n(i),size_n(i));
%     [~,~,~,cnt]=my_svd(A);
%     iter_time(i)=cnt;
% end
% figure;
% plot(size_n,iter_time);
% xlabel('size of n');
% ylabel('number of iterations');
% title('Number of iterations required by my\_svd for different size of n');
