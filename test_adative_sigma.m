% 测试收敛准则，自适应sigma
close all;clear all; clc
% a=[0 0;0 30;0,-30;30 0;-30 0;0 50;0 -50;50 0;-50 0];ks=9;%平移的坐标
a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
% a=[0 0;0 3;0,-3;3 0;-3 0;0 5;0 -5;5 0;-5 0];
% I0=fitsread('bbso_halph_fr_20150927_220002.fts');% KLL算法不会出现分歧
% I0=fitsread('bbso_halph_fl_20150414_173042.fts');% KLL算法不会出现分歧
I0=fitsread('bbso_halph_fl_20150924_210040.fts');% KLL算法会出现分歧
% I0=fitsread('bbso_halph_fl_20150924_204939.fts');% KLL算法会出现分歧
I0=double(I0(8:2032+7,8:2032+7));%size(a),figure;imshow(I0,[])
I0=I0-min(min(I0))+1;
F=double(imread('hd61900f.jpg'));
size(F),F0=F.*0.002+0.628;%figure;imshow(F,[]) %仿真平场，波幅范围0：0.512+0.628
% F0 = F0+0.01*rand(size(F,1));%对仿真下半场加1%随机噪声
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
    %         figure;imshow(D{1},[])
end
% I0=imcorr(D{1},D{6},15);
f=log(F0);
f= f(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
% figure;imshow(f,[])
% figure;imshow(D{1},[]);
O0=I0(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);




% MCC算法
tic
G0=D{1}.*0;
sig = 1.5;
K0=G0.*0;w=K0+1;n=w.*0+1e-10;
k=1;
for i=1:ks
    for j=i+1:ks
        if i==j
        else
            dij=-(a(i,:)-a(j,:));
            eijd=(D{i}-imshift(D{j},-dij(1),-dij(2))-G0+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
            qijd=exp(-eijd.^2./sig^2)+1e-10;%qijd=sig./qijd;
            % K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
            % K0=K0+qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
            K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2))...
                ;
            % n=n+qijd.*imshift(w,-dij(1),-dij(2));%+qijp.*imshift(w,dij(1),dij(2));
            % n=n+qijp.*imshift(w,dij(1),dij(2));%
            n=n+qijd.*imshift(w,-dij(1),-dij(2));
        end
    end
end
G0=K0./n;
first = std(G0(:));
TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%方差
sign(k) = sig ;
mccn(k) = TS;

while 1  %迭代次数
    k=k+1;
    K0=G0.*0;w=K0+1;n=w.*0+1e-10;
    for i=1:ks
        for j=i+1:ks
            if i==j
            else
                dij=-(a(i,:)-a(j,:));
                eijd=(D{i}-imshift(D{j},-dij(1),-dij(2))-G0+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
                qijd=exp(-eijd.^2./sig^2)+1e-10;%qijd=sig./qijd;
                % K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
                % K0=K0+qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
                K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2))...
                    ;
                % n=n+qijd.*imshift(w,-dij(1),-dij(2));%+qijp.*imshift(w,dij(1),dij(2));
                % n=n+qijp.*imshift(w,dij(1),dij(2));%
                n=n+qijd.*imshift(w,-dij(1),-dij(2));
            end
        end
    end
    G0=K0./n;
    scond = std(G0(:));
    if abs(scond-first)<10^-5
        break;
    end
    sig=sig+0.05*((scond-first)^2/sig^3-1/sig);
    TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%方差
    first=scond;
    sign(k) = sig ;
    mccn(k) = TS;
end
toc

figure,plot(sign)
figure,plot(mccn)

