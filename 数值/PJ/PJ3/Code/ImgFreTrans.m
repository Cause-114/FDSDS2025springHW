%% ImgFreTrans.m
% 实现图像在频域与时域之间的转化。
classdef ImgFreTrans

    methods (Access = public, Static)
% op代表是否对UV色度通道进行1：4采样。
% 如果op=1，那么输出上Y的长宽都将是UV的两倍。
% 将图像变化到频域上。
        function [Y,U,V] = img2fre(img,op)
            R = double(img(:, :, 1));
            G = double(img(:, :, 2));
            B = double(img(:, :, 3));
            Y = 0.299 * R + 0.587 * G + 0.114 * B;
            U = (B - Y)*0.564334; V = (R - Y)*0.713267;
            if(op==1)
                U=ImgFreTrans.sample(U);
                V=ImgFreTrans.sample(V);
            end
            Y = ImgFreTrans.t2f(Y-128);
            U = ImgFreTrans.t2f(U);
            V = ImgFreTrans.t2f(V);
        end
% 同上，这是将频域逆变换为图像的接口。
        function img = fre2img(Y,U,V,op)
            [m,n]=size(Y);
            Y = ImgFreTrans.f2t(Y)+128;
            U = ImgFreTrans.f2t(U);
            V = ImgFreTrans.f2t(V);
            if (op == 1)
                U = ImgFreTrans.dsample(U,m,n);
                V = ImgFreTrans.dsample(V,m,n);
            end
            B = U/0.564334 + Y; R = V/0.713267 + Y;
            G = (Y - 0.299 * R - 0.114 * B) / 0.587;
            img = uint8(round(cat(3, R, G, B)));
        end

    end

    methods (Access = private, Static)
% 将每个通道划分成一个个8*8通道后DCT变换到频域
        function X = t2f(x)
            [m, n] = size(x);
            X = zeros(m, n);
            for i = 1:8:m-7
                for j = 1:8:n-7
                    X(i:i + 7, j:j + 7) = DCT8by8.DCT_8by8(x(i:i + 7, j:j + 7));
                end
            end
        end
% 同上，这是逆变换。
        function x = f2t(X)
            [m, n] = size(X);
            x = zeros(m, n);
            for i = 1:8:m-7
                for j = 1:8:n-7
                    x(i:i + 7, j:j + 7) = DCT8by8.IDCT_8by8(X(i:i + 7, j:j + 7));
                end
            end
        end
% 对色度1：4采样，每4个点保留一个为4个点的平均值，这里加上了非偶数的边界条件处理
        function Us=sample(U)
            [m,n]=size(U);
            Us=zeros(round(m/2),round(n/2));
            for i=1:2:m
                for j=1:2:n
                    id=round(i/2);jd=round(j/2);
                    s=U(i,j);c=1;
                    if(i<m)
                        s=s+U(i+1,j);c=2;
                        if(j<m)
                            s=s+U(i,j+1)+U(i+1,j+1);
                            c=4;
                        end
                    elseif(j<m)
                        s=s+U(i,j+1);c=2;
                    end
                    Us(id,jd)=s/c;
                end
            end
        end
% 对采样后的还原。
        function U = dsample(Us,m,n)
            U = zeros(m, n);
            for i = 1:2:m
                for j = 1:2:n
                    id = round(i / 2); jd = round(j / 2);
                    U(i,j)=Us(id,jd);
                    if (i < m)
                        U(i+1,j)=Us(id,jd);
                        if (j < m)
                            U(i,j+1)=Us(id,jd);
                            U(i+1,j+1)=Us(id,jd);
                        end
                    elseif (j < m)
                        U(i,j+1)=Us(id,jd);
                    end
                end
            end
        end
    end
end
