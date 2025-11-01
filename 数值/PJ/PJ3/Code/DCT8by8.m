%% DCT8by8.m
% 提供8*8小块的DCT正变换与逆变换接口。
classdef DCT8by8
    properties (Access = private, Constant)
        A = DCT8by8.intial()
    end
    methods (Static)
        function A = intial()
            A = zeros(8, 8);
            for i = 2:8
                A(:, i) = cos(pi * (i - 1) / 8 * (1/2:1:15/2))';
            end
            A(:, 1) = 1 / sqrt(2);
            A = A ./ 2;
        end
        function x = DCT_8by8(x)
            x = DCT8by8.A' * x * DCT8by8.A;
        end
        function x = IDCT_8by8(x)
            x = DCT8by8.A * x * DCT8by8.A';
        end
    end
end
