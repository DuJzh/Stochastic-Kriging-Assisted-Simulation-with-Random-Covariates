function [ APFS ] = PFSEva3( Etestx,EGM,nd,delta,fi2,d,Bpred )
%PFSEva returns the APFS at a point Etestx 

xs = size(Etestx); tr = 0;

Emy = zeros(xs(1),nd); esvar = zeros(xs(1),nd);  
for temi = 1:nd
    Emy(:,temi) = SKpredictnew(EGM{temi},Etestx,Bpred);
    esvar(:,temi) = max(fi2{temi}(Etestx),0);  % matrix ill-conditioned
end
[bestval,bestdesig] = min(Emy,[],2);
nonval1 = Emy.'; nonvar1 = esvar.'; 
bestind = bestdesig + ((0:(xs(1)-1))*nd)';
bestvar = nonvar1(bestind);
nonval1(bestind) = [];
nonvar1(bestind) = [];
nonval = reshape(nonval1,(nd-1),xs(1)).';
nonvar = reshape(nonvar1,(nd-1),xs(1)).';
delta2 = nonval - bestval + delta;
inted = -delta2./sqrt(nonvar+bestvar);
APFS = sum(normcdf(inted),2);
if tr
    APFS = APFS';
end



end

