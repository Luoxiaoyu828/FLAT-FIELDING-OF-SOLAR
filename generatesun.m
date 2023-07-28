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