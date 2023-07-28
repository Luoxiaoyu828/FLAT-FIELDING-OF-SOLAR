close all;clear all; clc

% a=[0 0;0 30;0,-30;30 0;-30 0;0 50;0 -50;50 0;-50 0];ks=9;%平移的坐标
a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
% a=[0 0;0 3;0,-3;3 0;-3 0;0 5;0 -5;5 0;-5 0];

% I0=fitsread('bbso_halph_fr_20150927_220002.fts');%选择用于仿真的参考太阳图像
I0=fitsread('bbso_halph_fl_20150924_204939.fts');%选择用于仿真的参考太阳图像
I0=double(I0(8:2032+7,8:2032+7));%size(a),figure;imshow(I0,[])
I0=I0-min(min(I0))+1;
F=double(imread('hd61900f.jpg'));
size(F),F0=F.*0.002+0.628;%figure;imshow(F,[]) %仿真平场，波幅范围0：0.512+0.628
F0 = F0 + 0.01*rand(size(F,1));%对仿真下半场加1%随机噪声

figure,imshow(I0,[])
figure,imshow(F0,[])

[m,n]=size(F0);
for i=1:ks
    ai=a(i,:);
    F1=F0;%+0.01*rand(size(F,1));%对仿真下半场加1%随机噪声
    %   I{i}=generatesun(ai(1),ai(2),radius,sn).*F1; %生成系列太阳图像
    I1=imshift(I0,ai(1),ai(2));
    I1=I1.*F1;%+0.05*rand(size(F,1)).*max(max(I1));%对仿真图像加5%随机噪声
    I{i}=I1;%(round(m./2)-256:round(m./2)+256,round(n./2)-256:round(n./2)+256);%按指定图像生成序列图像
    %   figure;imshow(I{i},[])
end



radius=256;%模仿的太阳成像大小,取中间较均匀的区域进行平场估计
for i=1:ks
    I{i}= I{i}(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
    D{i}=log(I{i}); %取对数
    %         figure;imshow(D{i},[])
end
% I0=imcorr(D{1},D{6},15);
figure,imshow(I{2},[])
figure,imshow(I{4},[])
figure,imshow(I{7},[])
figure,imshow(I{10},[])
