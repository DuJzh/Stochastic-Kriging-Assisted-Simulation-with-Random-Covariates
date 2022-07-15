function [ y ] = Kernval( x,xcen,p,kernind )
%SqExpK returns the value of the squared exponential kernel
%   x is input in row xcen indicate the sampled covariate
% p is the parameter value
% t=0;
% if size(x,1) > 1
%     warning('x size abnormal')
%     x=x'; %t = 1;
% end
% 
% if size(x,1) > 1 && size(x,2) > 1
%     error('x size error')
% end
d1 = size(xcen,1); d2 = size(x,1); 

xi = repelem(x,d1,1); xc = repmat(xcen,d2,1);
disx = xi-xc; 
yt = sum(disx.^2,2);

distan = reshape(yt,d1,d2);  % each column corresponds to a x
% for di = 1:d2
%     distan(:,di) = sum( (x(di,:) - xcen).^2, 2 );
% end
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



% if t
%     y=y';
% end


end

