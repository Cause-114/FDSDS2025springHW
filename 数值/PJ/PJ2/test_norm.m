size_n = 50:250;
err = zeros(1, length(size_n));
iter_time=zeros(1,length(size_n));
for i = 1:length(size_n)
    A = rand(size_n(i), size_n(i));
    tic;
    n1=my_SSnorm2(A);
    iter_time(i)=toc;
    err(i) = (n1 - norm(A)) / norm(A);
end
figure;
plot(size_n, err);
xlabel('Size of matrix')
ylabel('Relative error')
title('Relative error of my\_SSnorm2()')

figure;
plot(log10(size_n),log10(iter_time));
xlabel('log10(size of n)');
ylabel('log10(time)');
title('Time required by my\_SSnorm2 for different matrix n*n');