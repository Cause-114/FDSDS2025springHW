%% Test.m
% 这是测试模块，提供对各个模块的主要功能的展示与评估的函数接口。
classdef Test
% 这四个是掩码矩阵，主要用于展示低频高频分量对图像的影响。
    properties (Access = private, Constant)
        Q0 = [1, 0, 0, 1, 0, 0, 1, 0; 1, 0, 0, 0, 1, 0, 0, 0; ...
                  0, 1, 1, 1, 0, 1, 0, 1; 0, 1, 0, 0, 1, 0, 1, 1; ...
                  1, 0, 1, 1, 1, 1, 0, 0; 0, 1, 1, 0, 1, 1, 1, 1; ...
                  0, 0, 1, 0, 0, 1, 0, 1; 0, 0, 0, 0, 0, 1, 0, 1; ];
        Q1 = [1, 1, 1, 0, 0, 0, 0, 0; 1, 1, 1, 0, 0, 0, 0, 0; ...
                  1, 1, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ];
        Q2 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 1, 1; ...
                  0, 0, 0, 0, 0, 1, 1, 1; 0, 0, 0, 0, 0, 1, 1, 1; ];
        Q3 = [1, 1, 1, 1, 1, 1, 1, 1; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ...
                  0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; ];
        Q4 = [0, 0, 1, 1, 1, 1, 1, 1; 0, 1, 1, 1, 1, 1, 1, 1; ...
                  1, 1, 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1, 1, 1; ...
                  1, 1, 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1, 1, 1; ...
                  1, 1, 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1, 1, 1; ];
    end

    methods (Access = public, Static)
% 给定一幅img:图片、plist:对应损失参数列表、op:是否进行UV通道1：4采样
% 输出对应损失参数下量化后频域YUV再霍夫曼编码下的所需要的总编码数（KB）。
        function percent = getCompressSiz(img, plist, op)
            percent = zeros(size(plist));
            [Y, U, V] = ImgFreTrans.img2fre(img, op);
            for i = 1:length(plist)
                [Y1, U1, V1] = FreFilter.Quantization(Y, U, V, plist(i));
                percent(i) = GetBit.getbit(Y1, U1, V1) / 8192;
            end
        end
% 参数类似于上一个函数，输出对应评估指标PSNR与SSIM
        function [psnr_l, ssim_l] = getQualityBand(img, ilist, op)
            psnr_l = zeros(size(ilist));
            ssim_l = zeros(size(ilist));
            [Y, U, V] = ImgFreTrans.img2fre(img, op);
            for i = 1:length(ilist)
                [Y1, U1, V1] = FreFilter.BandFilter(Y, U, V, 1, ilist(i));
                img1 = ImgFreTrans.fre2img(Y1, U1, V1, op);
                psnr_l(i) = psnr(img1, img);
                ssim_l(i) = ssim(img1, img);
            end
        end
% 参数类似于上一个函数，输出对应评估指标PSNR与SSIM
        function [psnr_l, ssim_l] = getQualityQuanti(img, plist, op)
            psnr_l = zeros(size(plist));
            ssim_l = zeros(size(plist));
            [Y, U, V] = ImgFreTrans.img2fre(img, op);
            for i = 1:length(plist)
                [Y1, U1, V1] = FreFilter.Quantization(Y, U, V, plist(i));
                [Y1, U1, V1] = FreFilter.DQuantization(Y1, U1, V1, plist(i));
                img1 = ImgFreTrans.fre2img(Y1, U1, V1, op);
                psnr_l(i) = psnr(img1, img);
                ssim_l(i) = ssim(img1, img);
            end
        end
% 参数类似于前面，展示量化-反量化操作对图片的影响
        function DisplayQuantization(img, op, p)
            [Y, U, V] = ImgFreTrans.img2fre(img, op);
            figure;
            if(nargin<3)
                for p1 = 1:4
                    [Y1, U1, V1] = FreFilter.Quantization(Y, U, V, p1);
                    [Y1, U1, V1] = FreFilter.DQuantization(Y1, U1, V1, p1);
                    img = ImgFreTrans.fre2img(Y1, U1, V1, op);
                    subplot(2, 2, p1); imshow(img);
                    title(sprintf("loss parameter p= %d", p1));
                end
            else
                [Y1, U1, V1] = FreFilter.Quantization(Y, U, V, p);
                [Y1, U1, V1] = FreFilter.DQuantization(Y1, U1, V1, p);
                img = ImgFreTrans.fre2img(Y1, U1, V1, op);
                imshow(img); title(sprintf("loss parameter p= %d", p));
            end
        end
% 参数类似于前面，展示高频、低频分量对图片的影响，依赖于前面的掩码矩阵
        function DisplayThreshold(img, op)
            [Y, U, V] = ImgFreTrans.img2fre(img, op);
            figure;
            subplot(3, 4, 1); imshow(ones(8, 8)); title('filter spectrum');
            subplot(3, 4, 2); imshow(img); title('result');
            [Y1, U1, V1] = FreFilter.MatrixFilter(Y, U, V, Test.Q0);
            img = ImgFreTrans.fre2img(Y1, U1, V1, op);
            subplot(3, 4, 3); imshow(Test.Q0);
            subplot(3, 4, 4); imshow(img);

            [Y1, U1, V1] = FreFilter.MatrixFilter(Y, U, V, Test.Q1);
            img = ImgFreTrans.fre2img(Y1, U1, V1, op);
            subplot(3, 4, 5); imshow(Test.Q1);
            subplot(3, 4, 6); imshow(img);

            [Y1, U1, V1] = FreFilter.MatrixFilter(Y, U, V, Test.Q2);
            img = ImgFreTrans.fre2img(Y1, U1, V1, op);
            subplot(3, 4, 7); imshow(Test.Q2);
            subplot(3, 4, 8); imshow(img);

            [Y1, U1, V1] = FreFilter.MatrixFilter(Y, U, V, Test.Q3);
            img = ImgFreTrans.fre2img(Y1, U1, V1, op);
            subplot(3, 4, 9); imshow(Test.Q3);
            subplot(3, 4, 10); imshow(img);

            [Y1, U1, V1] = FreFilter.MatrixFilter(Y, U, V, Test.Q4);
            img = ImgFreTrans.fre2img(Y1, U1, V1, op);
            subplot(3, 4, 11); imshow(Test.Q4);
            subplot(3, 4, 12); imshow(img);
        end

    end

end

