% This code is trying to show that my displayed DCT method which is based in 4N-FFT
% get the result exactly the same with the ordinary matrix-multiplication DCT.
N = 8;
a = rand(1, N);
X = DCT_FFT(a);
X1 = DCT(a);
fprintf("The norm of difference between DCT result by two different ways:%g\n",norm(X - X1));
fprintf("The norm of difference between rebulit DCT and intial sample:%g\n",norm(rebuilt(X) - a));

function x = DCT_FFT(a)
    % "Fast" DCT
    n = length(a);
    a1b = zeros(1, 4 * n);
    a1b(2:2:2 * n) = a;
    a1b = fft_cyy(a1b, 0);
    x = 2 / n * real(a1b(1:n));
    x(1) = x(1) / sqrt(2);
end

function X = DCT(x)
    % Matrix-multiplication based DCT, O(N^2)
    n = length(x);
    x = x(:).' * 2 / n;
    X = zeros(1, n);
    X(1) = sum(x) / sqrt(2);
    for i = 2:n
        X(i) = dot(x, cos(pi * (i - 1) / n * (1/2:1:n -1/2)));
    end
end

function y = rebuilt(x)
    n = length(x);
    y = ones(1, n)*x(1) / sqrt(2);
    for i = 2:n
        y = y + x(i) * cos((i - 1) * pi / n * (1/2:1:n -1/2));
    end
end

function x = fft_cyy(x, op)
    % Only for N=2^m length, i.e. N = 1,2,4...
    N = length(x);
    if (bitand(N, N - 1))
        error('N must be a power of 2');
    end
    x = x(:).';
    omega = exp(-2i * pi / N * (0:N / 2 - 1));
    if (op == 1)
        omega = conj(omega);
    end
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
