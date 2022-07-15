function [ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnS1( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu,Truva,Ttestx,Ttesty )
%TrnS returns the imse and ipfs under the truncated normal sampling ,Ttestx,Ttesty

nm = length(m); nd = size(desig,1); xd = size(desig,2);
% Tmu = span2+span1/2;
pd = makedist('Normal','mu',Tmu,'sigma',Truva); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);
TPmse1 = zeros(nm,1); TPpfs1 = zeros(nm,1);    % number of covariates layer
TPmse2 = zeros(nm,1); TPpfs2 = zeros(nm,1);    % number of covariates layer
TPmse3 = zeros(nm,1); TPpfs3 = zeros(nm,1);    % number of covariates layer
TPmse4 = zeros(nm,1); TPpfs4 = zeros(nm,1);    % number of covariates layer
for j = 1:nm    % number of covariates layer
    j
    Tmc = m(j); 
    TMmse1 = zeros(Lm,1); TMpfs1 = zeros(Lm,1);    % test covariate point layer
    TMmse2 = zeros(Lm,1); TMpfs2 = zeros(Lm,1);    % test covariate point layer
    TMmse3 = zeros(Lm,1); TMpfs3 = zeros(Lm,1);    % test covariate point layer
    TMmse4 = zeros(Lm,1); TMpfs4 = zeros(Lm,1);    % test covariate point layer

    for mi = 1:Lm    % number of random covariates layer
        tfun = fun;
        Tcontext = random(Tpd,Tmc,xd);      % random covarites generator
        Tms1 = zeros(nd,1); Tmy1 = zeros(testa,nd);  
        Tms2 = zeros(nd,1); Tmy2 = zeros(testa,nd);  
        Tms3 = zeros(nd,1); Tmy3 = zeros(testa,nd);  
        Tms4 = zeros(nd,1); Tmy4 = zeros(testa,nd);  
        
%         Ttestx = random(Tpd,testa,xd);         % test covariate generation for PFS
%         Ttesty = fun(desig,Ttestx);    % exact testfun value
        Targm = min(Ttesty,[],2);    % true optimal design
        
        for temi = 1:nd    % design layer
            Taccurval = fun(desig(temi,:),Tcontext); 
            TsampleCol = Taccurval+sd*randn(Tmc,n0);  %replications generator
            Tsample = mean(TsampleCol,2); tsd = sqrt(mean(var(TsampleCol,0,2))/n0);
            TGM1 = fitrgp(Tcontext,Tsample,'KernelFunction','exponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Tmy1(:,temi) = predict(TGM1,Ttestx);    % prediction of design temi at testx
            Tysd1 = (Tmy1(:,temi) - Ttesty(:,temi)).^2;
            Tms1(temi) = mean(Tysd1);
            
            TGM2 = fitrgp(Tcontext,Tsample,'KernelFunction','squaredexponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Tmy2(:,temi) = predict(TGM2,Ttestx);    % prediction of design temi at testx
            Tysd2 = (Tmy2(:,temi) - Ttesty(:,temi)).^2;
            Tms2(temi) = mean(Tysd2);
            
            TGM3 = fitrgp(Tcontext,Tsample,'KernelFunction','matern32','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Tmy3(:,temi) = predict(TGM3,Ttestx);    % prediction of design temi at testx
            Tysd3 = (Tmy3(:,temi) - Ttesty(:,temi)).^2;
            Tms3(temi) = mean(Tysd3);
            
            TGM4 = fitrgp(Tcontext,Tsample,'KernelFunction','matern52','FitMethod', ...
                'exact','PredictMethod','exact',...
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Tmy4(:,temi) = predict(TGM4,Ttestx);    % prediction of design temi at testx
            Tysd4 = (Tmy4(:,temi) - Ttesty(:,temi)).^2;
            Tms4(temi) = mean(Tysd4);
        end
        TMmse1(mi) = max(Tms1);
        TMmse2(mi) = max(Tms2);
        TMmse3(mi) = max(Tms3);
        TMmse4(mi) = max(Tms4);
        
        [~,Tapargmind1] = min(Tmy1,[],2);
        Tindtemp1 = (Tapargmind1-1)*testa + ((1:testa)'); 
        Tapargm1 = Ttesty(Tindtemp1);
        Tapres1 = abs(Tapargm1 - Targm) > Del;    %whether false selection happen at testx
        TMpfs1(mi) = mean(Tapres1);
        
        [~,Tapargmind2] = min(Tmy2,[],2);
        Tindtemp2 = (Tapargmind2-1)*testa + ((1:testa)'); 
        Tapargm2 = Ttesty(Tindtemp2);
        Tapres2 = abs(Tapargm2 - Targm) > Del;    %whether false selection happen at testx
        TMpfs2(mi) = mean(Tapres2);
        
        [~,Tapargmind3] = min(Tmy3,[],2);
        Tindtemp3 = (Tapargmind3-1)*testa + ((1:testa)'); 
        Tapargm3 = Ttesty(Tindtemp3);
        Tapres3 = abs(Tapargm3 - Targm) > Del;    %whether false selection happen at testx
        TMpfs3(mi) = mean(Tapres3);
        
        [~,Tapargmind4] = min(Tmy4,[],2);
        Tindtemp4 = (Tapargmind4-1)*testa + ((1:testa)'); 
        Tapargm4 = Ttesty(Tindtemp4);
        Tapres4 = abs(Tapargm4 - Targm) > Del;    %whether false selection happen at testx
        TMpfs4(mi) = mean(Tapres4);
    end
    TPmse1(j) = mean(TMmse1); TPpfs1(j) = mean(TMpfs1);    %average of random covariates level
    TPmse2(j) = mean(TMmse2); TPpfs2(j) = mean(TMpfs2);    %average of random covariates level
    TPmse3(j) = mean(TMmse3); TPpfs3(j) = mean(TMpfs3);    %average of random covariates level
    TPmse4(j) = mean(TMmse4); TPpfs4(j) = mean(TMpfs4);    %average of random covariates level
end

end

