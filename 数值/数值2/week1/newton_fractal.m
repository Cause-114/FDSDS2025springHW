lis_re = linspace(-2, 2, 500);
lis_im = linspace(-2, 2, 500);
colors=zeros(length(lis_re),length(lis_im));
for k=1:length(lis_re)
    for j=1:length(lis_im)
        z=lis_re(k)+lis_im(j)*1i;
        [t,z1]=get_iter_info(z);
        if real(z1)>0
            colors(j,k)=t-10;
        elseif imag(z1)>0
            colors(j,k)=t;
        else
            colors(j,k)=t+10;
        end
    end
end
imagesc(lis_re, lis_im, colors);
xlabel('Re(z)');
ylabel('Im(z)');
function [iter_time, z1] = get_iter_info(z)
    z1 = z - (z ^ 3 - 1) / (3 * z ^ 2);
    iter_time = -1;
    while (abs(z - z1) > 1e-6)
        iter_time = iter_time + 1;
        z = z1;
        z1 = z - (z ^ 3 - 1) / (3 * z ^ 2);
    end
end



% [x, y] = meshgrid(lis_re, lis_im);
% x = x(:); y = y(:);
% colors = zeros(length(x), 3);
% for k = 1:length(x)
%     z = x(k) + y(k) * 1i;
%     [t, z1] = get_iter_info(z);
%     if real(z1) > 0
%         colors(k, :) = [t / 9, 0, 0];
%     elseif imag(z1) > 0
%         colors(k, :) = [0, t / 9, 0];
%     else
%         colors(k, :) = [0, 0, t / 9];
%     end
% end
% figure;
% scatter(x, y, 2, colors, 'filled');