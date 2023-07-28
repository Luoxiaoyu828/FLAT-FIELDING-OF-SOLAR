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