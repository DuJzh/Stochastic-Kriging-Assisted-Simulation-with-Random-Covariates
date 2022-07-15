function [ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0 )
%UniS returns the imse and ipfs under the uniform sampling

nm = length(m); nd = size(desig,1); xd = size(desig,2); 
UPmse1 = zeros(nm,1); UPpfs1 = zeros(nm,1);    % number of covariates layer
UPmse2 = zeros(nm,1); UPpfs2 = zeros(nm,1);    % number of covariates layer
UPmse3 = zeros(nm,1); UPpfs3 = zeros(nm,1);    % number of covariates layer
UPmse4 = zeros(nm,1); UPpfs4 = zeros(nm,1);    % number of covariates layer
for j = 1:nm    % number of covariates layer
    j
    Umc = m(j);
    UMmse1 = zeros(Lm,1); UMpfs1 = zeros(Lm,1);    % test covariate point layer
    UMmse2 = zeros(Lm,1); UMpfs2 = zeros(Lm,1);    % test covariate point layer
    UMmse3 = zeros(Lm,1); UMpfs3 = zeros(Lm,1);    % test covariate point layer
    UMmse4 = zeros(Lm,1); UMpfs4 = zeros(Lm,1);    % test covariate point layer

    for mi = 1:Lm    % number of random covariates layer
        tfun = fun;
        Ucontext = rand(Umc,xd)*span1+span2;    % random covarites generator
        Ums1 = zeros(nd,1); Umy1 = zeros(testa,nd);  
        Ums2 = zeros(nd,1); Umy2 = zeros(testa,nd);  
        Ums3 = zeros(nd,1); Umy3 = zeros(testa,nd);  
        Ums4 = zeros(nd,1); Umy4 = zeros(testa,nd);  
        
        Utestx = rand(testa,xd)*span1+span2;   % test covariate generation for PFS
        Utesty = fun(desig,Utestx);    % exact testfun value
        Uargm = min(Utesty,[],2);    % true optimal design
        
        for temi = 1:nd    % design layer
            Uaccurval = fun(desig(temi,:),Ucontext); 
            UsampleCol = Uaccurval+sd*randn(Umc,n0);  %replications generator
            Usample = mean(UsampleCol,2); tsd = sqrt(mean(var(UsampleCol,0,2))/n0);
            UGM1 = fitrgp(Ucontext,Usample,'KernelFunction','exponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Umy1(:,temi) = predict(UGM1,Utestx);    % prediction of design temi at testx
            Uysd1 = (Umy1(:,temi) - Utesty(:,temi)).^2;
            Ums1(temi) = mean(Uysd1);
            
            UGM2 = fitrgp(Ucontext,Usample,'KernelFunction','squaredexponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Umy2(:,temi) = predict(UGM2,Utestx);    % prediction of design temi at testx
            Uysd2 = (Umy2(:,temi) - Utesty(:,temi)).^2;
            Ums2(temi) = mean(Uysd2);
            
            UGM3 = fitrgp(Ucontext,Usample,'KernelFunction','matern32','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Umy3(:,temi) = predict(UGM3,Utestx);    % prediction of design temi at testx
            Uysd3 = (Umy3(:,temi) - Utesty(:,temi)).^2;
            Ums3(temi) = mean(Uysd3);
            
            UGM4 = fitrgp(Ucontext,Usample,'KernelFunction','matern52','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Umy4(:,temi) = predict(UGM4,Utestx);    % prediction of design temi at testx
            Uysd4 = (Umy4(:,temi) - Utesty(:,temi)).^2;
            Ums4(temi) = mean(Uysd4);
        end
        UMmse1(mi) = max(Ums1);
        UMmse2(mi) = max(Ums2);
        UMmse3(mi) = max(Ums3);
        UMmse4(mi) = max(Ums4);
        
        [~,Uapargmind1] = min(Umy1,[],2);
        Uindtemp1 = (Uapargmind1-1)*testa + ((1:testa)'); 
        Uapargm1 = Utesty(Uindtemp1);
        Uapres1 = abs(Uapargm1 - Uargm) > Del;    %whether false selection happen at testx
        UMpfs1(mi) = mean(Uapres1);
        
        [~,Uapargmind2] = min(Umy2,[],2);
        Uindtemp2 = (Uapargmind2-1)*testa + ((1:testa)'); 
        Uapargm2 = Utesty(Uindtemp2);
        Uapres2 = abs(Uapargm2 - Uargm) > Del;    %whether false selection happen at testx
        UMpfs2(mi) = mean(Uapres2);
        
        [~,Uapargmind3] = min(Umy3,[],2);
        Uindtemp3 = (Uapargmind3-1)*testa + ((1:testa)'); 
        Uapargm3 = Utesty(Uindtemp3);
        Uapres3 = abs(Uapargm3 - Uargm) > Del;    %whether false selection happen at testx
        UMpfs3(mi) = mean(Uapres3);
        
        [~,Uapargmind4] = min(Umy4,[],2);
        Uindtemp4 = (Uapargmind4-1)*testa + ((1:testa)'); 
        Uapargm4 = Utesty(Uindtemp4);
        Uapres4 = abs(Uapargm4 - Uargm) > Del;    %whether false selection happen at testx
        UMpfs4(mi) = mean(Uapres4);
    end
    UPmse1(j) = mean(UMmse1); UPpfs1(j) = mean(UMpfs1);    %average of random covariates level
    UPmse2(j) = mean(UMmse2); UPpfs2(j) = mean(UMpfs2);    %average of random covariates level
    UPmse3(j) = mean(UMmse3); UPpfs3(j) = mean(UMpfs3);    %average of random covariates level
    UPmse4(j) = mean(UMmse4); UPpfs4(j) = mean(UMpfs4);    %average of random covariates level
end

end

