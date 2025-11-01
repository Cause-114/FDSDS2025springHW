lis = linspace(-50, 50, 1000);
res0 = sin(lis); % actual ans.
res1 = zeros(1, length(lis)); % store the result of method (a).
res2 = zeros(1, length(lis)); % store the result of method (b).

for i = 1:length(lis)
    x = lis(i);
    res1(i) = comp_sin_seri(1e-6, x);
    xx = Shift2Norm(x);
    res2(i) = comp_sin_seri(1e-6, xx);
    % you can change the error require, I set 1e-6 above.
end
figure;
plot(lis, log10(abs(res0-res1)));
ylabel("log10 error")
title("using the method (a)-log10 version of truncation error")
figure;
plot(lis, log10(abs(res0-res2)));
ylabel("log10 error")
title("using the method (b)-log10 version of truncation error")
% figure
% plot(lis,res1);
% figure
% plot(lis,res2);
function xx = Shift2Norm(x)
    xx = mod(x, 2 * pi);
    if xx > 1.5 * pi
        xx = xx - 2 * pi;
    elseif xx > 0.5 * pi
        xx = pi - xx;
    end
end

function res = comp_sin_seri(err, x)
    res = 0; ser = x; i = 3;
    while (abs(ser) > err)
        res = res + ser;
        ser = ser * x ^ 2 * -1 / (i * (i - 1));
        i = i + 2;
    end
end
