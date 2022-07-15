function [ y ] = Kernval( x,xcen,p,kernind )
%Kernval returns the value of the kernel
% x is the covariates that we want to calculate the value of kernel
% xcen is the covariate points already sampled
% p is the parameter value
% kernind=1: exponential; kernind=2: squared-exponential; 
% kernind=3: matern 3/2; kernind=4: matern 5/2; 

d1 = size(xcen,1); d2 = size(x,1); 

xi = repelem(x,d1,1); xc = repmat(xcen,d2,1);
disx = xi-xc; 
yt = sum(disx.^2,2);

distan = reshape(yt,d1,d2);  % each column corresponds to a x

p2 = p.^2;

if kernind == 1  %exponential
    y = p2(2)*exp(-sqrt(distan)/p(1));
elseif kernind == 2  % squared-exponential
    y = p2(2)*exp(-distan/p2(1)/2);
elseif kernind == 3  % matern 3/2
    rs = sqrt(3)*sqrt(distan)/p(1);
    y = p2(2)*(1+rs).*exp(-rs);
elseif kernind == 4  % matern 5/2
    rs = sqrt(5)*sqrt(distan)/p(1);
    y = p2(2)*(1+rs+rs.^2/3).*exp(-rs);
end



end

