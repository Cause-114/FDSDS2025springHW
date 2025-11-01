%% GetBit.m
% 这个类用以得到理论上的JEPG位流总位数（bit）
classdef GetBit

    properties (Access = private, Constant)
% ACL1,ACL2,ACC1,ACC2两个一组对应AC分量的亮度、色度行程对（n,L）霍夫曼树简化形式
% 因为其实在n，L较大时，霍夫曼编码长度几乎都是16，所以可以将n分为较小，较大两种情形，从而减少一大片的16.
% 以色度（ACC1，ACC2）为例，在n>3,L>3时，行程对霍夫曼编码长度全是16，因此只需要将n>3的数据只存储L=1,L=2的情形为ACC2。
% 而将n=0,1,2,3的所有L全存为ACC1（同时将EOB,ZRL两种情况编入(0,0)与(1,0),将所有数据往后移一位）。
        ACL1 = [4, 2, 2, 3, 4, 5, 7, 8, 10, 16, 16; ...
                    11, 4, 5, 7, 9, 11, 16, 16, 16, 16, 16];
        ACL2 = [5, 8, 10, 12; 6, 9, 12, 16; 6, 10, 16, 16; 7, 11, 16, 16; ...
                    7, 12, 16, 16; 8, 12, 16, 16; 9, 15, 16, 16; 9, 16, 16, 16; ...
                    9, 16, 16, 16; 10, 16, 16, 16; 10, 16, 16, 16; 11, 16, 16, 16; ...
                    16, 16, 16, 16; 16, 16, 16, 16];
        ACC1 = [2, 2, 3, 4, 5, 5, 6, 7, 9, 10, 12; 10, 4, 6, 8, 9, 11, 12, 16, 16, 16, 16; ...
                    0, 5, 8, 10, 12, 15, 16, 16, 16, 16, 16; 0, 5, 8, 10, 12, 16, 16, 16, 16, 16, 16];
        ACC2 = [6, 9; 6, 10; 7, 11; 7, 11; 8, 16; 9, 16;
                9, 16; 9, 16; 9, 16; 11, 16; 14, 16; 15, 16];
% 8*8小块之字形行动顺序。
        zigzag = [1, 2, 6, 7, 15, 16, 28, 29, ...
                      3, 5, 8, 14, 17, 27, 30, 43, ...
                      4, 9, 13, 18, 26, 31, 42, 44, ...
                      10, 12, 19, 25, 32, 41, 45, 54, ...
                      11, 20, 24, 33, 40, 46, 53, 55, ...
                      21, 23, 34, 39, 47, 52, 56, 61, ...
                      22, 35, 38, 48, 51, 57, 60, 62, ...
                      36, 37, 49, 50, 58, 59, 63, 64];
    end

    methods (Access = public, Static)
% 这个是对外开放的接口，从Y,U,V三个通道获取霍夫曼编码对应的总字节数。
        function res = getbit(Y, U, V)
            res = 0;
        % 亮度Y
            [m, n] = size(Y);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    res = res + GetBit.proL(Y(i:i + 7, j:j + 7));
                end
            end
        % 色度UV
            [m, n] = size(U);
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    res = res + GetBit.proC(U(i:i + 7, j:j + 7));
                end
            end
            for i = 1:8:m - 7
                for j = 1:8:n - 7
                    res = res + GetBit.proC(V(i:i + 7, j:j + 7));
                end
            end
        end

    end

    methods (Access = private, Static)
% 处理亮度
        function res = proL(X)
            X = round(X(GetBit.zigzag));
            n = 0; L = 0;
            if (X(1) ~= 0)
                L = floor(log2(abs(X(1)))) + 1;
            end
            res = GetBit.DCLHuffLen(L); tmp = 0;
            % 按之字形行动
            for i = 2:64
                if (X(i) == 0)
                    n = n + 1;
                    if (n == 15)
                        tmp = tmp + GetBit.ACLHuffLen(1, 0); n = 0;
                    % 连续遇到15个0，这个时候需要引入ZRL中断，但是为例避免是EOB的情况，我们将其暂计入tmp.
                    elseif (i == 64)
                        res = res + GetBit.ACLHuffLen(0, 0);
                    % EOB.即以一串0结尾。
                    end
                else
                    % 说明还未结尾，将之前还未记录的ZRL的编码总和以及这次的编码长度加入res。
                    L = floor(log2(abs(X(i)))) + 1;
                    res = res + tmp + GetBit.ACLHuffLen(n, L);
                    n = 0; tmp = 0;
                end
            end
        end
% 处理色度
        function res = proC(X)
            X = round(X(GetBit.zigzag));
            n = 0; L = 0;
            if (X(1) ~= 0)
                L = floor(log2(abs(X(1)))) + 1;
            end
            res = GetBit.DCCHuffLen(L); tmp = 0;
            for i = 2:64
                if (X(i) == 0)
                    n = n + 1;
                    if (n == 15)
                        tmp = tmp + GetBit.ACCHuffLen(1, 0); n = 0;
                    elseif (i == 64)
                        res = res + GetBit.ACCHuffLen(0, 0);
                    end
                else
                    L = floor(log2(abs(X(i)))) + 1;
                    res = res + tmp + GetBit.ACCHuffLen(n, L);
                    n = 0; tmp = 0;
                end
            end
        end
% 返回亮度DC分量的某一尺寸所需编码长度（前缀霍夫曼编码长度+实际数据所需长度L）
        function res = DCLHuffLen(L)
% 查表能够发现有这样的规律
            if (L == 0)
                res = 2;
            elseif (L > 4)
                res = 2 * L - 2;
            else
                res = 3 + L;
            end
        end
% 同上，这是对色度处理（霍夫曼表不一样）
        function res = DCCHuffLen(L)
% 查表能够发现有这样的规律
            if (L < 2)
                res = 2 + L;
            else
                res = L * 2;
            end
        end
% 返回亮度AC分量的某一尺寸所需编码长度（前缀霍夫曼编码长度+实际数据所需长度L）
        function res = ACLHuffLen(n, L)
% 具体逻辑请见对properties的说明。 
            if (n <= 1)
                res = GetBit.ACL1(n + 1, min(L + 1, 11));
            else
                if (L > 4)
                    res = 16;
                else
                    res = GetBit.ACL2(min(n - 1, 14), L);
                end
            end
            res = res + L;
        end
% 返色度AC分量的某一尺寸所需编码长度（表不一样，所以重写）
        function res = ACCHuffLen(n, L)
            if (n < 4)
                res = GetBit.ACC1(n + 1, min(L + 1, 11));
            else
                if (L > 2)
                    res = 16;
                else
                    res = GetBit.ACC2(min(n - 3, 12), L);
                end
            end
            res = res + L;
        end
    end

end