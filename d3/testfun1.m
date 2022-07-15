function [ y ] = testfun1( x,xcen )
%testfun return the mean value of De Jong's function
%   x is input; xcen indicate the covariate

d = size(xcen,1);
dx = size(x,1);


xi = repelem(x,d,1); xc = repmat(xcen,dx,1);
disx = xi-xc; 
yt = sum(disx.^2,2);
y = reshape(yt,[d,dx]);


end

