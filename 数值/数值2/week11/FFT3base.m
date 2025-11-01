test_fft_cyy3(1e-10, 1000, 3.^(0:5));
figure_maker();

function x = fft_cyy3(x)
    % Only for N=3^m length, i.e. N = 1,3,9...
    N = length(x);
    if (3 ^ round(log2(N) / log2(3)) ~= N)
        error('N must be a power of 3');
    end
    x = x(:).';
    omega = exp(-2i * pi / N * (0:2 * N / 3));
    step = N / 3;
    while (step >= 1)
        for j = 1:step
            series1 = j:3 * step:N; series2 = j + step:3 * step:N;
            series3 = j + 2 * step:3 * step:N;
            x(series2) = x(series2) .* omega(1:step:N / 3);
            x(series3) = x(series3) .* omega(1:2 * step:2 * N / 3);
            tmp1 = x(series1) + x(series2) + x(series3);
            tmp2 = x(series1) + omega(N / 3 + 1) * x(series2) + omega(2 * N / 3 + 1) * x(series3);
            tmp3 = x(series1) + omega(2 * N / 3 + 1) * x(series2) + omega(N / 3 + 1) * x(series3);
            x(j:step:N) = [tmp1, tmp2, tmp3];
        end
        step = step / 3;
    end
end

function test_fft_cyy3(tol, test_cases, Len)
    for j = 1:length(Len)
        flag = 1;
        for k = 1:test_cases
            x = rand(1, Len(j));
            n = norm(fft_cyy3(x) - fft(x));
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

function figure_maker()
    l = 3 .^ (1:12);
    time_cost = zeros(1, length(l));
    for i = 1:length(l)
        x = rand(1, l(i));
        tic;
        fft_cyy3(x);
        time_cost(i) = log2(toc / i)/log2(3);
    end
    figure;
    plot(1:12, time_cost, 'o-');
    xlabel('log3 version of input length');
    ylabel('log3(Time cost(s)/log3(N))');
    title('relation between time cost and N of our fft3 function');
end