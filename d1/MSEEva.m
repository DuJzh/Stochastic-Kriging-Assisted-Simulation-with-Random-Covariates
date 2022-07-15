


function [ y ] = MSEEva( testx,Econshort,KernParac,coinvc,FSFc,FSc,fbasis,kernind,nd )
%MSEEva returns the MSE 

xs = size(testx); mse = zeros(xs(1),nd);

for temi = 1:nd
    KernPara = KernParac{temi};
    msep1 = KernPara(2)^2; 
    coinv = coinvc{temi};
    fi1 = Kernval(testx,Econshort,KernPara,kernind);
    etav = fbasis{temi}(testx) - FSc{temi}*fi1;
    mse(:,temi) = msep1-sum(fi1'*coinv .*fi1', 2)...
        +sum(etav'*FSFc{temi}.* etav', 2);
end
y = max(mse,[],2);

end





