function [ output_args ] = new_kll_huai( input_args )
close all;clear all

% % a=[0 0;349   155; -368  127;-130 150;-130 -150;-200 -200;100 -200;250 -200;150 00];ks=9;%ƽ�Ƶ�����
% % a=[0 0;0 30;0,-30;30 0;-30 0;0 50;0 -50;50 0;-50 0];ks=9;%ƽ�Ƶ�����
% a=[0,0;50,-14; 40,17; 21,43; -5,56;-28,35;-44,7;-51,-24;-27,40;5,-43;37,-36];ks=11;
% % a=[0 0;0 3;0,-3;3 0;-3 0;0 5;0 -5;5 0;-5 0];ks=9
% % a=[0 0;1 0;0,-1;1 0;-1 0;1 1;-1 1;1 -1;-1 -1];ks=9
% 
% % I0=fitsread('bbso_halph_fr_20150927_220002.fts');%ѡ�����ڷ���Ĳο�̫��ͼ��
% I0=fitsread('bbso_halph_fl_20150924_204939.fts');%ѡ�����ڷ���Ĳο�̫��ͼ��
% I0=double(I0(8:2032+7,8:2032+7));%size(a),figure;imshow(I0,[])
% I0=I0-min(min(I0))+1;
% F=double(imread('hd61900f.jpg'));
% size(F),F0=F.*0.002+0.628;%figure;imshow(F,[]) %����ƽ����������Χ0��0.512+0.628
% [m,n]=size(F0);
% for i=1:ks
%   ai=a(i,:);
%   F1=F0;%+0.01*rand(size(F,1));%�Է����°볡��1%�������
% %   I{i}=generatesun(ai(1),ai(2),radius,sn).*F1; %����ϵ��̫��ͼ��
%   I1=imshift(I0,ai(1),ai(2));
%   I1=I1.*F1+0.05*rand(size(F,1)).*max(max(I1));%�Է���ͼ���5%�������
%   I{i}=I1;%(round(m./2)-256:round(m./2)+256,round(n./2)-256:round(n./2)+256);%��ָ��ͼ����������ͼ��
% %   figure;imshow(I{i},[])
% end
% 
% % figure;imshow(I{1},[]);
% % figure;imshow(I{3},[]);
% % a,
% % 
clear all
sc=1;%������ʵ��������ݽ���ƽ������
I{1}=fitsread('HrHa061229053108full.fit','Primary');
I{2}=fitsread('HrHa061229053129flat.fit','Primary');
I{3}=fitsread('HrHa061229053148flat.fit','Primary');
I{4}=fitsread('HrHa061229053159flat.fit','Primary');
I{5}=fitsread('HrHa061229053208flat.fit','Primary');
I{6}=fitsread('HrHa061229053229flat.fit','Primary');
I{7}=fitsread('HrHa061229053238flat.fit','Primary');
I{8}=fitsread('HrHa061229053318flat.fit','Primary');
I{9}=fitsread('HrHa061229053348flat.fit','Primary');
I{10}=fitsread('HrHa061229053108full.fit','Primary');
% I{1}=fitsread('HrHa110606032059FULL.fit','Primary');
% I{2}=fitsread('HrHa110606032338FULL.fit','Primary');
% I{3}=fitsread('HrHa110606032323FULL.fit','Primary');
% I{4}=fitsread('HrHa110606032305FULL.fit','Primary');
% I{5}=fitsread('HrHa110606032244FULL.fit','Primary');
% I{6}=fitsread('HrHa110606032226FULL.fit','Primary');
% I{7}=fitsread('HrHa110606032208FULL.fit','Primary');
% I{8}=fitsread('HrHa110606032153FULL.fit','Primary');
% I{9}=fitsread('HrHa110606032129FULL.fit','Primary');
% I{10}=fitsread('HrHa110606032059FULL.fit','Primary');
figure;imshow((I{10}),[]);%return
for i=1:length(I)
        I{i} = swapbytes(int16(I{i})) ;
        I{i} = double(I{i}(:,1:end-8)) ;
%         figure,imshow(I{i},[]) ;
end
% figure;imshow(log(I{10}),[]);
[m,n]=size(I{1}),
for i=1:9   %�Զ�������������ƽ�����꣬�Է���ͼ����Ҫ
    cc = fftshift(ifft2(fft2(I{1}).* fft2(I{i})));
%     imshow(cc,[])
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
a(i,:)=[xpeak, ypeak];
if i>1,a(i,:)=[(a(1,1)-xpeak), (a(1,2)-ypeak)];end
end
c=a(1,1),yc=a(1,2);
a(1,:)=[0,0];
ks=9;F0=I{10};
I0=log(I{1});f=log(I{10});


for i=1:ks    
wb{i}=I{i}>70;  %ȷ��ʵ��ͼ��ı߽�
% figure;imshow(I{i}>50,[]);
end
% a,

tic
itnumber=15;%��������
radius=500;%ģ�µ�̫�������С,ȡ�м�Ͼ��ȵ��������ƽ������
for i=1:ks
    I{i}= I{i}(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
    D{i}=log(I{i}+2); %����ȡ������ȡ����
%         figure;imshow(D{i},[])
end

f=log(F0);
f= f(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
figure;imshow(f,[])  %��ʾ��׼��ƽ��
figure;imshow(D{1},[]); %��ʾ���ĵ�ͼ��
O0=I0(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);
for i=1:9
wb{i}=wb{i}(round(m./2)-radius:round(m./2)+radius,round(n./2)-radius:round(n./2)+radius);    
end

%Kuhn method����׼KUHN�㷨
G0=D{1}.*0;K0=G0;
% w=G0+1;
n=G0;
for i=1:ks
    for j=1:ks
        if i==j
        else
        dij=-(a(i,:)-a(j,:));
%         K0=K0+D{i}-imshift(D{j},-dij(1),-dij(2))+D{j}-imshift(D{i},dij(1),dij(2));%δ���Ǳ߽�Ӱ��
        K0=K0+(D{i}-imshift(D{j},-dij(1),-dij(2))).*wb{i}.*imshift(wb{j},-dij(1),-dij(2))...
            +(D{j}-imshift(D{i},dij(1),dij(2))).*wb{j}.*imshift(wb{i},dij(1),dij(2));%���Ǳ߽�Ӱ��
%         n=n+imshift(w,-dij(1),-dij(2))+imshift(w,dij(1),dij(2));%δ���Ǳ߽�Ӱ��
        n=n+wb{i}.*imshift(wb{j},-dij(1),-dij(2))+wb{j}.*imshift(wb{i},dij(1),dij(2));        
        end
    end
end
K0=K0./n;
G1=K0;
% figure;imshow(G1,[]);
G1(isnan(G1))=0;w=G1~=0;%figure;imshow(w)
% return
for k=1:itnumber  %iterative G(r+1)=K+sum(G(r))��������
G0=K0.*0;n=G0;
for i=1:ks
    for j=i+1:ks
       if i==j
        else
        dij=-(a(i,:)-a(j,:));
        G0=G0+imshift(G1,-dij(1),-dij(2)).*imshift(w,-dij(1),-dij(2))+imshift(G1,dij(1),dij(2)).*imshift(w,dij(1),dij(2));        
%         G0=G0+imshift(G1,-dij(1),-dij(2))+imshift(G1,dij(1),dij(2));  
        n=n+imshift(w,-dij(1),-dij(2))+imshift(w,dij(1),dij(2));
       end
    end
end
G1=G0./n+K0;
G1(isnan(G1))=0;w=G1~=0;
% figure;imshow(G1,[])
end
figure;imshow(G1,[])   
% figure;plot(G1(radius,:));hold on;plot(f(radius,:),'r');
TT=G1-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), TS=sqrt(mean(TTS(:))),std(abs(TT(:)))
% return

% %JONGCHUL CHAE 2004 method
% F=D{1}.*0;
% O1=F;w=F+1;
% for i=1:ks
%     dxy=(a(i,:)-a(1,:));
%     O1=O1+imshift(D{i},-dxy(1),-dxy(2));
% end
% O=O1./ks;
% [m,n]=size(O);
% for i=1:ks
% C(i)=sum(D{i}(:))./(m.*n) -sum(O(:))./(m.*n);
% end
% % figure;imshow(O,[])
% % C
% % return
% 
% for i=1:itnumber
% on=w.*0+1e-10;fn=on;Df=O.*0;Do=Df;Dk=Do;
% for k=1:ks
%     dxy=-(a(k,:)-a(1,:));
% %     dxy=-a(k,:);
%     Dk=(imshift(O,-dxy(1),-dxy(2))+F-D{k}+C(k)).*imshift(w,-dxy(1),-dxy(2));
%     Df=Df+(imshift(O,-dxy(1),-dxy(2))+F-D{k}+C(k)).*imshift(w,-dxy(1),-dxy(2));
%     fn=fn+imshift(w,-dxy(1),-dxy(2));
%     Do=Do+(C(k)+O+imshift(F,dxy(1),dxy(2))-imshift(D{k},dxy(1),dxy(2))).*imshift(w,dxy(1),dxy(2));
%     on=on+imshift(w,dxy(1),dxy(2));
% %     C(k)=C(k)-sum(Dk(:))./sum(sum(imshift(w,-dxy(1),-dxy(2))));
% end
% O=O-Do./on;
% F=F-Df./fn;
% % C,
% end
% F=F-sum(F(:))./(m.*n);
% OF=O-sum(O(:))./(m.*n);
% % figure;plot(F(256,:));hold on;plot(F0(256,:),'r');
% figure;imshow(F,[]);figure;imshow(OF,[]);
% TT=F-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), TS=sqrt(mean(TTS(:))),std(TT(:))%����
% TT=O-log(O0);TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), %JC method 2004 
% 
% %MCC-JONGCHUL CHAE 2004 method new method
% F=D{1}.*0;
% O1=F;w=F+1;
% for i=1:ks
%     dxy=(a(i,:)-a(1,:));
%     O1=O1+imshift(D{i},-dxy(1),-dxy(2));
% end
% O=O1./ks;
% [m,n]=size(O);
% for i=1:ks
% C(i)=sum(D{i}(:))./(m.*n) -sum(O(:))./(m.*n);
% end
% figure;imshow(O,[])
% % C
% % return
% for i=1:itnumber
% on=w.*0+1e-10;fn=on;Df=O.*0;Do=Df;Dk=Do;sig=0.005;
% for k=1:ks
%     dxy=-(a(k,:)-a(1,:));
%     eijf=(imshift(O,-dxy(1),-dxy(2))+F-D{k}).*imshift(w,-dxy(1),-dxy(2));
%     qijf=100*exp(-eijf.^2./sig)+1e-10;
%     eijo=(imshift(F,dxy(1),dxy(2))+O-imshift(D{k},dxy(1),dxy(2))).*imshift(w,dxy(1),dxy(2));
%     qijo=100*exp(-eijo.^2./sig)+1e-10;
% 
%     Df=Df+qijf.*(D{k}-imshift(O,-dxy(1),-dxy(2))).*imshift(w,-dxy(1),-dxy(2));
%     fn=fn+qijf.*imshift(w,-dxy(1),-dxy(2));
%     Do=Do+qijo.*(imshift(D{k},dxy(1),dxy(2))-imshift(F,-dxy(1),-dxy(2))).*imshift(w,dxy(1),dxy(2));
%     on=on+qijo.*imshift(w,dxy(1),dxy(2));
% end
% O=Do./on;
% F=Df./fn;
% % C,
% end
% F=F-sum(F(:))./(m.*n);
% OF=O-sum(O(:))./(m.*n);
% % figure;plot(F(256,:));hold on;plot(F0(256,:),'r');
% figure;imshow(F,[]);figure;imshow(OF,[]);
% TT=F-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), TS=sqrt(mean(TTS(:))),std(TT(:))%����
% TT=O-log(O0);TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), %JC method 2004 
% % return

% MCC-based Kuhn method ��KUHN�㷨
% figure;plot(F(256,:))
G0=D{1}.*0;w=G0+1;
for k=1:itnumber  %��������
n=w.*0+1e-10;K0=G0.*0;sig=0.005;
for i=1:ks
    for j=i+1:ks
        if i==j
        else
        dij=-(a(i,:)-a(j,:));
        eijd=(D{i}-imshift(D{j},-dij(1),-dij(2))-G0+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
        qijd=100*exp(-eijd.^2./sig)+1e-10;%qijd=sig./qijd;
        eijp=imshift(D{i},dij(1),dij(2))-D{j}-imshift(G0,dij(1),dij(2))+G0;
        qijp=100*exp(-eijp.^2./sig)+1e-10;%qijp=sig./qijp;
%         K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2));
%         K0=K0+qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
        K0=K0+qijd.*(D{i}-imshift(D{j},-dij(1),-dij(2))+imshift(G0,-dij(1),-dij(2))).*imshift(w,-dij(1),-dij(2))...
        +qijp.*(-imshift(D{i},dij(1),dij(2))+D{j}+imshift(G0,dij(1),dij(2))).*imshift(w,dij(1),dij(2));
%         n=n+qijd.*imshift(w,-dij(1),-dij(2));%+qijp.*imshift(w,dij(1),dij(2));
%         n=n+qijp.*imshift(w,dij(1),dij(2));%
        n=n+qijd.*imshift(w,-dij(1),-dij(2))+qijp.*imshift(w,dij(1),dij(2));
        end
    end
end
G0=K0./n;
% if k>=1
% figure;plot(K0(256,:));%imshow(K0,[])
% end
end
G0(isnan(G0))=0;
figure;imshow(G0,[]);
figure;plot(G0(radius,:));%MCC-Kuhn method
hold on;plot(f(radius,:),'r');%standard flatfielding
hold on;plot(G1(radius,:),'g');%standard kuhn method
% hold on;plot(F(radius,:),'r');%JC method 2004
TT=G0-f;TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), TS=sqrt(mean(TTS(:))),std(abs(TT(:)))%����

% figure;imshow(O,[]); %Log object of JC method 2004 
figure;imshow(D{1}-G0,[]);%Log object of MCC-Kuhn method 
figure;imshow(D{1}-G1,[]);%Log object of standard kuhn method 
figure;imshow(log(O0),[]);%Log object of original data 
% TT=O-log(O0);TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))), %JC method 2004 
TT=D{1}-G0-log(O0);TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))),%MCC-Kuhn method 
TT=D{1}-G1-log(O0);TTS=TT.^2;psnr=10*log10(100./mean(TTS(:))),%standard kuhn method
toc
figure,imshow(exp(G0),[]),title('MCCƽ��')
figure,imshow(exp(G1),[]),title('KLLƽ��')
figure,imshow(exp(D{1}-G0),[]),title('MCCĿ��')
figure,imshow(exp(D{1}-G1),[]),title('KLLĿ��')
figure,imshow(exp(O0),[]),title('ԭʼĿ��')
return


function I=generatesun(x0,y0,radius,sn)
%generage a sun image with center at (x0,y0)
%the image size is about [1024,1024],and the radius 
for x=-256:256
    for y=-256:256
        r=sqrt((x+x0).^2+(y+y0).^2);sita=asin(r./radius);
        I(x+257,y+257)=0;
        if r<=radius && sita<=pi./2
            I(x+257,y+257)=10000*(1-0.88-0.23+0.88.*cos(sita)+0.23*cos(sita).^2);
            if I(x+257,y+257)<0; I(x+257,y+257)=0;end
        end
    end
end
I=I+100;%add 50 as the sky background
I1=I.*0;I1(200:400,300:500)=I1(200:400,300:500)+100;I=I+I1;%add 
% m=max(I(:)); I = imnoise(double(I)./m,'gaussian',0,sn/10000).*m;
I=I+sn*randn(size(I,1));%add random noise with amplitude sn
return

function I0=imshift(I0,x,y)
%the center of image I0 is shifted from (0,0) to (x,y),x,y are integers
%Left/right 
s = size(I0);
shiftp = x; 
I0 = I0(:,mod((1:s(2))+s(2)+shiftp,s(2))+1);
%Up/Down 
shiftp = y; 
I0 = I0(mod((1:s(1))+s(1)+shiftp,s(1))+1,:); 
return


function I0=imcorr(x,y,w)
%compute the correlation of X,Y with window w
%Left/right 
[m,n]=size(x);
I0=x.*0;
for i=w+1:m-w
    for j=w+1:n-w
        a=x(i-w+1:i+w);b=y(i-w+1:i+w);
        c=a.*b;
        I0(i,j)=sum(c(:))/sqrt(sum(dot(a,a))*sum(dot(b,b)));
    end
end
I0=I0';
return