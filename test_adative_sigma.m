% ��������׼������Ӧsigma
close all;clear all; clc
% a=[0 0;0 30;0,-30;30 0;-30 0;0 50;0 -50;50 0;-50 0];ks=9;%ƽ�Ƶ�����
a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
% a=[0 0;0 3;0,-3;3 0;-3 0;0 5;0 -5;5 0;-5 0];
% I0=fitsread('bbso_halph_fr_20150927_220002.fts');% KLL�㷨������ַ���
% I0=fitsread('bbso_halph_fl_20150414_173042.fts');% KLL�㷨������ַ���
I0=fitsread('bbso_halph_fl_20150924_210040.fts');% KLL�㷨����ַ���
% I0=fitsread('bbso_halph_fl_20150924_204939.fts');% KLL�㷨����ַ���
I0=double(I0(8:2032+7,8:2032+7));%size(a),figure;imshow(I0,[])
I0=I0-min(min(I0))+1;
F=double(imread('hd61900f.jpg'));
size(F),F0=F.*0.002+0.628;%figure;imshow(F,[]) %����ƽ����������Χ0��0.512+0.628
% F0 = F0+0.01*rand(size(F,1));%�Է����°볡��1%�������
[m,n]=size(F0);
for i=1:ks
    ai=a(i,:);
    F1=F0;%+0.01*rand(size(F,1));%�Է����°볡��1%�������
    %   I{i}=generatesun(ai(1),ai(2),radius,sn).*F1; %����ϵ��̫��ͼ��
    I1=imshift(I0,ai(1),ai(2));
    I1=I1.*F1;%+0.05*rand(size(F,1)).*max(max(I1));%�Է���ͼ���5%�������
    I{i}=I1;%(round(m./2)-256:round(m./2)+256,round(n./2)-256:round(n./2)+256);%��ָ��ͼ����������ͼ��
    %   figure;imshow(I{i},[])
end



radius=256;%ģ�µ�̫�������С,ȡ�м�Ͼ��ȵ��������ƽ������
for i=1:ks
    I{i}= I{i}(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
    D{i}=log(I{i}); %ȡ����
    %         figure;imshow(D{1},[])
end
% I0=imcorr(D{1},D{6},15);
f=log(F0);
f= f(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
% figure;imshow(f,[])
% figure;imshow(D{1},[]);
O0=I0(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);




% MCC�㷨
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
TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%����
sign(k) = sig ;
mccn(k) = TS;

while 1  %��������
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
    TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%����
    first=scond;
    sign(k) = sig ;
    mccn(k) = TS;
end
toc

figure,plot(sign)
figure,plot(mccn)

