% ÕºœÒ≈˙¡øœ‘ æ
close all;clear all; clc
str1 = dir(['D:\TDDownload\tongji\*.fts']);
for xu_kk =1:length(str1)
    clearvars -except str1 xu_kk M
    temp = str1(xu_kk).name;
    temp=temp(21:22);
    riqi{xu_kk}=temp;
    filenamee = fullfile('D:\TDDownload\tongji',str1(xu_kk).name);
    a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
    I0=fitsread(filenamee);
    if size(I0,1)~=size(I0,2)
        disp('the rows of raw image are not the same with the cols')
    end
    I0 = double(imresize(I0,[256,256]));
    figure,imshow(I0,[])
    title(num2str(xu_kk));
    M(xu_kk) = getframe;
end
% movie(M)