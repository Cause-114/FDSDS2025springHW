%% main.m

% 值得说明的是，正如实验报告中所提到的那样，由于我们是将RGB转为YUV处理的，并且人眼对色度UV敏感度低于
% Y通道的亮度。所以实际主流压缩处理偏向于对UV通道1：4采样再压缩。你能够看到这个main文件调用的几乎所有
% 函数带一个“1”的参数。这代表接受对UV通道的1：4采样。如果你不想要这样，请将对应的函数的这个参数改成0.

% 你也可以在这里换成其他你喜欢的图片，但是最好保持像素点能整除8，因为我们处理实际上依赖8*8小块。
% 注意路径，这里读取文件是相对路径，即从当前目录下的TestFigure文件夹读取图片。我交上去的版本有两个测试PNG。
img = imread("./TestFigure/Lenna.png");

% 展示两种方式的压缩效果，第一种是用掩码矩阵将零位置的频域值给抹去。第二种是量化操作。这里也可以单独显示
% 某一个损失参数下的复原图片，使用方法见注释掉的第15行(在第三个参数传入你想要的p)。
Test.DisplayThreshold(img, 1);
Test.DisplayQuantization(img, 1);
% Test.DisplayQuantization(img, 1, 1);

% 得到在不同压缩参数下对图像质量的评估参数。
% 第一种是带状滤波即将频域上横纵坐标之和大于某项的给置0，这个主要是用作对比，用以
% 说明PSNR、SSIM确实是越大图像质量越高。
% 第二种就是不同损失参数、有无对色度通道采样的下的量化-反量化得到的图像。
% [psnr_l1, ssim_l1] = Test.getQualityBand(img, 1:15, 1);
[psnr_l1, ssim_l1] = Test.getQualityQuanti(img, 1:8, 0);
[psnr_l2, ssim_l2] = Test.getQualityQuanti(img, 1:8, 1);

% 得到理论上的基于霍夫曼编码的JPEG位流总位数（已经除了8192，直接是对应的KB）
siz = Test.getCompressSiz(img, 1:8, 1);
siz1 = Test.getCompressSiz(img, 1:8, 0);
% 以512*512*3/1024=768KB作为未压缩的参考大小。
percent = 1 - (siz ./ (numel(img) / 1024));

% 画出对量化的两种评估指标与损失参数之间的关系。
figure; plot(1:8, psnr_l1, 1:8, psnr_l2);
legend("No sample of UV", "1:4 sample of UV");
title("PSNR VALUE of Quantization")
ylabel("PNSR(dB)"); xlabel("loss parameter p");
figure; plot(1:8, ssim_l1, 1:8, ssim_l2);
legend("No sample of UV", "1:4 sample of UV");
title("SSIM VALUE of Quantization")
ylabel("SSIM"); xlabel("loss parameter p");

% % 画出对带状滤波（主要是对比）的两种评估指标与损失参数之间的关系。
% figure; plot(1:15, psnr_l1);
% title("PSNR VALUE of Band Filter")
% ylabel("PNSR(dB)"); xlabel("loss parameter p");
% figure; plot(1:15, ssim_l1);
% title("SSIM VALUE of Band Filter")
% ylabel("SSIM"); xlabel("loss parameter p");

% 画出JPEG总位流大小（KB）与损失参数之间的关系。
figure;
plot(1:8, siz1, 1:8, siz);
legend("No sample of UV", "1:4 sample of UV");
ylabel("Storage space(KB)"); xlabel("loss parameter");
title("theoretical Storage space with loss parameter")

% 画出大小与质量随损失参数的图像（同一张图下）
figure;
yyaxis left;
plot(1:8, percent);
ylabel("Reduced Storage percent by Quantization");
xlabel("loss parameter p");
yyaxis right;
plot(1:8, ssim_l2);
title("Evaluate of Compress (storage vs effect)")
ylabel("SSIM");
