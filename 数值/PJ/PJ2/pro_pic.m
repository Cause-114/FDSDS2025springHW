img=imread("pictures\zby.jpg");
disp(size(img));
rchannel=double(img(:,:,1));
gchannel=double(img(:,:,2));
bchannel=double(img(:,:,3));
[Ur,Vr,Sr]=my_svd(rchannel);
[Ug,Vg,Sg]=my_svd(gchannel);
[Ub,Vb,Sb]=my_svd(bchannel);
rank_list=[18,44,88,131,175];
for i = 1:length(rank_list)
    r=rank_list(i);
    name=['zby_r' num2str(r) '.jpg'];
    img_low_r=uint8(Ur(:,1:r)*Sr(1:r,1:r)*Vr(:,1:r)');
    img_low_g=uint8(Ug(:,1:r)*Sg(1:r,1:r)*Vg(:,1:r)');
    img_low_b=uint8(Ub(:,1:r)*Sb(1:r,1:r)*Vb(:,1:r)');
    img_low=cat(3,img_low_r,img_low_g,img_low_b);
    imwrite(img_low,['pictures\' name]);
end
