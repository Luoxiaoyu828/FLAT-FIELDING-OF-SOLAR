% 统计结果
close all;clear all; clc
str1 = dir(['D:\TDDownload\tongji\*.fts']);
for xu_kk =1:length(str1)
    clearvars -except str1 xu_kk jc mcc2 kll riqi time_mcc time_kll time_jc
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
    % figure,imshow(I0,[])
    I0=I0-min(min(I0))+1;
    F=imread('hd61900f.jpg');
    if size(F,1)~=size(F,2)
        disp('the rows of flat are not the same with the cols')
    end
    F = double(imresize(F,[256,256]));
    F0=F.*0.002+0.628;%figure;imshow(F,[]) %仿真平场，波幅范围0：0.512+0.628
    % F0 = F0+0.01*rand(size(F,1));%对仿真下半场加1%随机噪声
    [m,n]=size(F0);
    for i=1:ks
        ai=a(i,:);
        F1=F0;%+0.01*rand(size(F,1));%对仿真下半场加1%随机噪声
        %   I{i}=generatesun(ai(1),ai(2),radius,sn).*F1; %生成系列太阳图像
        I1=imshift(I0,ai(1),ai(2));
        I1=I1.*F1;%+0.05*rand(size(F,1)).*max(max(I1));%对仿真图像加5%随机噪声
        I{i}=I1;%(round(m./2)-256:round(m./2)+256,round(n./2)-256:round(n./2)+256);%按指定图像生成序列图像
        D{i}=log(I{i}); %取对数
    end
    f=log(F0);
    % figure;imshow(f,[])
    % figure;imshow(D{1},[]);
    O0=I0;
    
    % MCC算法
    tic;
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
        sig=sig+0.05*((scond-first)^2/sig^3-1/sig);
        TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));%std(TT(:))%方差
        first=scond;
        sign(k) = sig ;
        mccn(k) = TS;
        if abs(scond-first)<10^-5
            break;
        end
    end
    time_2=toc;
    mcc2(xu_kk) = TS;
    time_mcc(xu_kk) = time_2;
    %     figure,plot(sign)
    %     figure,plot(mccn)
    figure,subplot(1,2,1),imshow(G0,[]),title('mcc')
    subplot(1,2,2),imshow(f,[]),title('raw image')
    
    % Kuhn算法
    % 先生成K(x)
    itnumber=15;
    tic;
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
%                 temp_2 = G1-f;TTS=temp_2.^2; TS=sqrt(mean(TTS(:)));
%                 temp_xu(k) = TS;
    end
    time_2=toc;
    time_kll(xu_kk) = time_2;
%         figure,plot(temp_xu)
    TT=G1-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:)));
    kll(xu_kk) = TS ;
    
    % JC方法
    tic;
    itnumber=15;
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
%                 temp_1 = F-sum(F(:))./(m.*n);
%                 temp_2 = temp_1-f;TTS=temp_2.^2; TS=sqrt(mean(TTS(:)));
%                 temp_xu(i) = TS;
    end
%         figure,plot(temp_xu)
    F=F-sum(F(:))./(m.*n);
    OF=O-sum(O(:))./(m.*n);
    TT=F-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))); TS=sqrt(mean(TTS(:))); %方差
    jc(xu_kk) = TS ;
    time_2=toc;
    time_jc(xu_kk) = time_2;
end

figure,plot(jc,'r')
hold on,plot(kll,'k')
hold on,plot(mcc2,'b')
legend('JC','KLL','MCC')

figure,plot(time_jc,'r')
hold on,plot(time_kll,'k')
hold on,plot(time_mcc,'b')
legend('JC','KLL','MCC')

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');
for i=1:length(riqi)
    iter(i)=str2num(riqi{i});
end
xlim([0,31])
ylim([0.2,0.33])
X1 = [iter', iter', iter'] ;
YMatrix1 = [kll', jc', mcc2'] ;
loglog1 = plot(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','+','LineStyle','--','DisplayName','Kuhn et al.');
set(loglog1(2),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','*','LineStyle','--','DisplayName','Chae J');
set(loglog1(3),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','o','DisplayName','Zheng et al.');
xlabel('Date','FontWeight','bold','FontName','Times New Roman');
ylabel('Standard Error(\epsilon)','FontWeight','bold','FontName','Times New Roman');
legend(axes1,'show');

% time
figure1 = figure;
axes1 = axes('Parent',figure1,...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');
for i=1:length(riqi)
    iter(i)=str2num(riqi{i});
end
xlim([0,31])
ylim([0.5,5.5])
X1 = [iter', iter', iter'] ;
YMatrix1 = [time_kll', time_jc', time_mcc'] ;
loglog1 = plot(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','+','LineStyle','--','DisplayName','Kuhn et al.');
set(loglog1(2),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','*','LineStyle','--','DisplayName','Chae J');
set(loglog1(3),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','o','DisplayName','Zheng et al.');
xlabel('Date','FontWeight','bold','FontName','Times New Roman');
ylabel('Time(s)','FontWeight','bold','FontName','Times New Roman');
legend(axes1,'show');