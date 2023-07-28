function disp = algin(image1,image2)
% image = imread('h619l.jpg') ;
% image = imread('h619l03.jpg') ;
% figure,imshow(image)
% bw = im2bw(image,0.2) ;
% figure,imshow(bw,[])
% bw2 = edge(bw,'sobel') ;
% bw2(:,1:10) = 0 ;
% bw2(:,end-10:end) = 0 ;
% bw2(1:10,:) = 0 ;
% bw2(end-10:end,:) = 0 ;
% [x,y] = find(bw2==1) ;
% figure,imshow(bw2,[])
% hold on,plot(y,x,'.r')
% x=x(:);
% y=y(:);
% m=[x y ones(size(x))]\[-(x.^2+y.^2)];
% xc = -.5*m(1)%拟合圆心X轴数值
% yc = -.5*m(2)%拟合圆心Y轴数值
% R  =  sqrt((m(1)^2+m(2)^2)/4-m(3))%拟合半径数值
% xc,yc,R
% figure,imshow(image) ;
% viscircles([yc,xc], R,'EdgeColor','b');
[xc1,yc1] = circlefit(image1) ;
[xc2,yc2] = circlefit(image2) ;
disp = [yc1-yc2,xc1-xc2] ;
end

function [xc,yc] = circlefit(image)
image = uint8(image) ;
bw = im2bw(image,0.2) ;
bw2 = edge(bw,'sobel') ;
bw2(:,1:10) = 0 ;
bw2(:,end-10:end) = 0 ;
bw2(1:10,:) = 0 ;
bw2(end-10:end,:) = 0 ;
[x,y] = find(bw2==1) ;
x=x(:);
y=y(:);
m=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*m(1); % 拟合圆心X轴数值
yc = -.5*m(2); % 拟合圆心Y轴数值
R  =  sqrt((m(1)^2+m(2)^2)/4-m(3)); % 拟合半径数值
% figure,imshow(image) ;
% viscircles([yc,xc], R,'EdgeColor','b');
end
