% 收敛情况对比
close all;clear all; clc
% a=[0 0;0 30;0,-30;30 0;-30 0;0 50;0 -50;50 0;-50 0];ks=9;%平移的坐标
a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
% a=[0 0;0 3;0,-3;3 0;-3 0;0 5;0 -5;5 0;-5 0];
I0=fitsread('bbso_halph_fr_20150927_220002.fts');% KLL算法不会出现分歧
% I0=fitsread('bbso_halph_fl_20150414_173042.fts');% KLL算法不会出现分歧
% I0=fitsread('bbso_halph_fl_20150924_210040.fts');% KLL算法会出现分歧
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
figure;imshow(f,[])
figure;imshow(D{1},[]);
O0=I0(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);

% itnumber=15;%迭代次数
xut = [1:19, 20:5:70] ;
errn = zeros(size(xut));
for xu_itnumber=1:length(xut)
    itnumber = xut(xu_itnumber) 
    % Kuhn算法
    % 先生成K(x)
    tic
    G0=D{1}.*0;K0=G0;
    w=G0+1;n=G0;
    for i=1:ks
        for j=1:ks
            if i==j
            else
                dij=-(a(i,:)-a(j,:));
                K0=K0+D{i}-imshift(D{j},-dij(1),-dij(2))+D{j}-imshift(D{i},dij(1),dij(2));
                n=n+imshift(w,-dij(1),-dij(2))+imshift(w,dij(1),dij(2));
            end
        end
    end
    K0=K0./n;
    G1=0;
    for k=1:itnumber  %iterative G(r+1)=K+sum(G(r))迭代次数
        G0=K0.*0;n=G0+1e-10;
        for i=1:ks
            for j=i+1:ks
                if i==j
                else
                    dij=-(a(i,:)-a(j,:));
                    G0=G0+imshift(G1,-dij(1),-dij(2)).*imshift(w,-dij(1),-dij(2))+imshift(G1,dij(1),dij(2)).*imshift(w,dij(1),dij(2));
                    n=n+imshift(w,-dij(1),-dij(2))+imshift(w,dij(1),dij(2));
                end
            end
        end
        G1=G0./n+K0;
    end
    toc
    TT=G1-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), TS=sqrt(mean(TTS(:)));
    kll(itnumber) = TS ;

    % JC方法
    F=D{1}.*0;
    O1=F;w=F+1;
    for i=1:ks
        dxy=(a(i,:)-a(1,:));
        % O1=O1+imshift(D{i},-dxy(1),-dxy(2));
        O1=O1+ D{i};
    end
    O=O1./ks;
    [m,n]=size(O);
    for i=1:ks
        C(i)=sum(D{i}(:))./(m.*n) -sum(O(:))./(m.*n);
    end
    for i=1:itnumber
        on=w.*0+1e-10;fn=on;Df=O.*0;Do=Df;Dk=Do;
        for k=1:ks
            dxy=-(a(k,:)-a(1,:));
            %     dxy=-a(k,:);
            Dk=(imshift(O,-dxy(1),-dxy(2))+F-D{k}+C(k)).*imshift(w,-dxy(1),-dxy(2));
            Df=Df+(imshift(O,-dxy(1),-dxy(2))+F-D{k}+C(k)).*imshift(w,-dxy(1),-dxy(2));
            fn=fn+imshift(w,-dxy(1),-dxy(2));
            Do=Do+(C(k)+O+imshift(F,dxy(1),dxy(2))-imshift(D{k},dxy(1),dxy(2))).*imshift(w,dxy(1),dxy(2));
            on=on+imshift(w,dxy(1),dxy(2));
            C(k)=C(k)-sum(Dk(:))./sum(sum(imshift(w,-dxy(1),-dxy(2))));
        end
        O=O-Do./on;
        F=F-Df./fn;
    end
    F=F-sum(F(:))./(m.*n);
    OF=O-sum(O(:))./(m.*n);
    TT=F-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:))); %方差
    jc(itnumber) = TS ;

    % MCC算法
    tic    
    G0=D{1}.*0;
    % sign = 1.5:-0.01:0.01;
    sign = 0.8:-0.01:0.01;
    mccn = zeros(size(sign));    
    % for k=1:length(sign)  %迭代次数
    for k=1:itnumber  %迭代次数
        sig = sign(k);
        n=w.*0+1e-10;K0=G0.*0;w=K0+1;
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
        temp22=K0./n;
        n=w.*0+1e-10;K0=G0.*0;
        for i=1:ks
            for j=i+1:ks
                if i==j
                else
                    dij=-(a(i,:)-a(j,:));
                    eijp=imshift(D{i},dij(1),dij(2))-D{j}-imshift(G0,dij(1),dij(2))+G0;
                    qijp=exp(-eijp.^2./sig^2)+1e-10;%qijp=sig./qijp;
                    % K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
                    % K0=K0+qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
                    K0=K0+qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
                    % n=n+qijd.*imshift(w,-dij(1),-dij(2));%+qijp.*imshift(w,dij(1),dij(2));
                    % n=n+qijp.*imshift(w,dij(1),dij(2));%
                    n=n+qijp.*imshift(w,dij(1),dij(2));
                end
            end
        end
        temp11=K0./n;
        G0 = (temp22+temp11)/2 ;
        TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%方差
        mccn(k) = TS ;
    end
    toc
    mccn(mccn==0) = 5000;
    [TS,t_index] = min(mccn);%方差
    errn(itnumber) = sign(t_index);
    mcc2(itnumber) = TS ;

end
% figure,plot(1:length(kll),kll,'-')
% hold on,plot(1:length(mcc2),mcc2,'--r')
iter = find(kll>0) ;
figure,plot(kll(iter),'-r')
hold on,plot(jc(iter),'-b')
hold on,plot(mcc2(iter),'-b')
figure,plot(errn(iter));
title('sigma');
save('jingdu.mat','kll', 'mcc2','jc')