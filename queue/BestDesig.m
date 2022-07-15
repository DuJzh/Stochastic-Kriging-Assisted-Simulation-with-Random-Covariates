function [ y ] = BestDesig( x,xcen,cu,U )
%BestDesig return the mean value of each setting (lambda, x).
%   x is design; xcen indicate the covariate
% cu: service cost
% U: Upper bound of total cost

d = size(xcen,1);
dx = size(x,1);

xi = repelem(x,d,1); xc = repmat(xcen,dx,1);
yt1 = 1./(xi-xc) + cu .* xi; yt = min(yt1, U);
y = reshape(yt,[d,dx]);


end

