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