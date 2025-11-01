f = @(x) arrayfun(@(x) square_wave(x), x);
figure_maker(f, 2, 2 .^ (0:5));

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

function py = trInter(y, p)
    % It's quite noticable that as I use my self-defined fft_cyy function, the
    % interpolation requires length(y) to be a power of 2.
    % In general, if you replace fft_cyy with fft, you can use any length(y) that is odd.
    N = length(y);
    y = fft_cyy(y, 0);
    y = y / N;
    py = zeros(1, p);
    py(1:floor(N / 2 + 1)) = y(1:floor(N / 2 + 1));
    py(p - floor(N / 2) + 2:p) = y(floor(N / 2 + 2):N);
    py = real(fft_cyy(py, 1));
    fprintf("p = %d, N = %d\n", p, N);
end

function y = square_wave(x)

    if (mod(x, 1) < 1e-10)
        y = 0;
    elseif (mod(x, 2) < 1)
        y = 1;
    else
        y = -1;
    end

end

function figure_maker(f, T, N)
    % From left to right: fucntion handel of periodic signal,
    % period T, and length of samples, i.e. N=[N1,N2,N3].
    px = linspace(0, T, 512);
    l = length(N);
    figure;
    plot(px, f(px)); hold on;
    legends = cell(1, l + 1);
    legends{1} = "Original signal";

    for i = 1:l
        y = f((0:N(i) - 1) * T / N(i));
        py = trInter(y, 512);
        plot(px, py); hold on;
        legends{i + 1} = sprintf("Interpolated signal, N = %d", N(i));
    end

    title("Interpolation of a square wave");
    xlabel("Time");
    ylabel("value");
    legend(legends);
end
