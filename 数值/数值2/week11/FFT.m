test_fft_cyy(1e-12, 1000, 2.^(0:5));
figure_maker();

function x = fft_cyy(x)
% Only for N=2^m length, i.e. N = 1,2,4...
    N = length(x);
    if (bitand(N, N - 1))
        error('N must be a power of 2');
    end
    x = x(:).';
    omega = exp(-2i * pi / N * (0:N / 2 - 1));
    step = bitshift(N, -1);
    while (step >= 1)
        for j = 1:step
            x(j + step:2 * step:N) = x(j + step:2 * step:N) .* omega(1:step:N / 2);
            tmp = x(j:2 * step:N) + x(j + step:2 * step:N);
            x(j + N / 2:step:N) = x(j:2 * step:N) - x(j + step:2 * step:N);
            x(j:step:j + N / 2 - 1) = tmp;
        end
        step = bitshift(step, -1);
    end
end

function figure_maker()
    l = 2 .^ (1:16);
    time_cost = zeros(1, length(l));
    for i = 1:length(l)
        x = rand(1, l(i));
        tic;
        fft_cyy(x);
        time_cost(i) = log2(toc / i);
    end
    figure;
    plot(1:16, time_cost, 'o-');
    xlabel('log2 version of input length');
    ylabel('log2(Time cost(s)/log2(N))');
    title('relation between time cost and N of our fft function');
end

function test_fft_cyy(tol, test_cases, Len)
    for j = 1:length(Len)
        flag = 1;
        for k = 1:test_cases
            x = rand(1, Len(j));
            n = norm(fft_cyy(x) - fft(x));
            if (n > tol)
                fprintf('The difference between the two results is:%g\n', n);
                flag = 0;
                break;
            end
        end
        if (flag)
            fprintf("Under tolerance:%g, %d %d-length ", tol, test_cases, Len(j))
            fprintf("test cases, our function is same as matlab's fft function.\n");
        end
    end
end
