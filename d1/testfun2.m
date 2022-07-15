function [ y ] = testfun2( x,xcen )
%testfun return the mean value of Griewangk's function
%   x is input; xcen indicate the covariate

[d,L] = size(xcen); 
dx = size(x,1);


xi = repelem(x,d,1); xc = repmat(xcen,dx,1);
disx = xi-xc; 
dif = repmat(sqrt(1:L), d*dx , 1); csdisx = cos( disx./dif );
yt = 1/4000*sum(disx.^2,2) - prod(csdisx,2) + 1;
y = reshape(yt,[d,dx]);


end

