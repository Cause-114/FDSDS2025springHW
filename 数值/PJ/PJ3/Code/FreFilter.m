%% FreFilter.m
% 对频域上的YUV通道进行简单压缩。
classdef FreFilter
% 这两个是JPEG标准中推荐的亮度、色度量化矩阵的取值。
    properties (Access = private, Constant)
        Qy = [16, 11, 10, 16, 24, 40, 51, 61; 12, 12, 14, 19, 26, 58, 60, 55; ...
                  14, 13, 16, 24, 40, 57, 69, 56; 14, 17, 22, 29, 51, 87, 80, 62; ...
                  18, 22, 37, 56, 68, 109, 103, 77; 24, 35, 55, 64, 81, 104, 113, 92; ...
                  49, 64, 78, 87, 103, 121, 120, 101; 72, 92, 95, 98, 112, 100, 103, 99];
        Qc = [17, 18, 24, 47, 99, 99, 99, 99; 18, 21, 26, 66, 99, 99, 99, 99; ...
                  24, 26, 56, 99, 99, 99, 99, 99; 47, 66, 99, 99, 99, 99, 99, 99; ...
                  99, 99, 99, 99, 99, 99, 99, 99; 99, 99, 99, 99, 99, 99, 99, 99; ...
                  99, 99, 99, 99, 99, 99, 99, 99; 99, 99, 99, 99, 99, 99, 99, 99];
    end

    methods (Access = public, Static)
% 带状滤波。具体来讲，保留横纵坐标之和加一在s到e之间的频域值。
        function [Y, U, V] = BandFilter(Y, U, V, s, e)
            Y = FreFilter.bandfilter(Y, s, e);
            U = FreFilter.bandfilter(U, s, e);
            V = FreFilter.bandfilter(V, s, e);
        end
% 利用掩码矩阵二值化频域。
        function [Y, U, V] = MatrixFilter(Y, U, V, Q)
            Y = FreFilter.matrixfilter(Y, Q);
            U = FreFilter.matrixfilter(U, Q);
            V = FreFilter.matrixfilter(V, Q);
        end
% 量化
        function [Y, U, V] = Quantization(Y, U, V, p)
            Y = FreFilter.quantization(Y, FreFilter.Qy * p);
            U = FreFilter.quantization(U, FreFilter.Qc * p);
            V = FreFilter.quantization(V, FreFilter.Qc * p);
        end
% 反量化
        function [Y, U, V] = DQuantization(Y, U, V, p)
            Y = FreFilter.dquantization(Y, FreFilter.Qy * p);
            U = FreFilter.dquantization(U, FreFilter.Qc * p);
            V = FreFilter.dquantization(V, FreFilter.Qc * p);
        end

    end
% 这下面是对每个通道划分成小块分而治之，不对外开放接口。
    methods (Access = private, Static)

        function X = bandfilter(X, s, e)
            [m, n] = size(X);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    for it = 0:7
                        s1 = min(max(s - it - 2, -1), 7);
                        e1 = max(min(e - it, 8), 0);
                        X(i + it, j + [0:s1, e1:7]) = 0;
                    end
                end
            end
        end

        function X = matrixfilter(X, Q)
            [m, n] = size(X);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    X(i:i + 7, j:j + 7) = X(i:i + 7, j:j + 7) .* Q;
                end
            end
        end

        function X = quantization(X, Q)
            [m, n] = size(X);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    X(i:i + 7, j:j + 7) = round(X(i:i + 7, j:j + 7) ./ Q);
                end
            end
        end

        function X = dquantization(X, Q)
            [m, n] = size(X);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    X(i:i + 7, j:j + 7) = X(i:i + 7, j:j + 7) .* Q;
                end
            end
        end

    end

end
