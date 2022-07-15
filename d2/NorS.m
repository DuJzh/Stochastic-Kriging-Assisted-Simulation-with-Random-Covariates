function [ NPmse1,NPpfs1,NPmse2,NPpfs2,NPmse3,NPpfs3,NPmse4,NPpfs4 ] ...
    = NorS( fun,m,desig,scen,ssd,Lm,testa,Del,sd,n0,ssd2 )
%TrnS returns the imse and ipfs under the truncated normal sampling

nm = length(m); nd = size(desig,1); xd = size(desig,2);
NPmse1 = zeros(nm,1); NPpfs1 = zeros(nm,1);    % number of covariates layer
NPmse2 = zeros(nm,1); NPpfs2 = zeros(nm,1);    % number of covariates layer
NPmse3 = zeros(nm,1); NPpfs3 = zeros(nm,1);    % number of covariates layer
NPmse4 = zeros(nm,1); NPpfs4 = zeros(nm,1);    % number of covariates layer
for j = 1:nm    % number of covariates layer
    j
    Nmc = m(j); 
    NMmse1 = zeros(Lm,1); NMpfs1 = zeros(Lm,1);    % test covariate point layer
    NMmse2 = zeros(Lm,1); NMpfs2 = zeros(Lm,1);    % test covariate point layer
    NMmse3 = zeros(Lm,1); NMpfs3 = zeros(Lm,1);    % test covariate point layer
    NMmse4 = zeros(Lm,1); NMpfs4 = zeros(Lm,1);    % test covariate point layer

    for mi = 1:Lm    % number of random covariates layer
        tfun = fun;
        Ncontext = randn(Nmc,xd)*ssd2+scen;   % random covarites generator
        Nms1 = zeros(nd,1); Nmy1 = zeros(testa,nd);  
        Nms2 = zeros(nd,1); Nmy2 = zeros(testa,nd);  
        Nms3 = zeros(nd,1); Nmy3 = zeros(testa,nd);  
        Nms4 = zeros(nd,1); Nmy4 = zeros(testa,nd);  
        
        Ntestx = randn(testa,xd)*ssd+scen;         % test covariate generation for PFS
        Ntesty = fun(desig,Ntestx);    % exact testfun value
        Nargm = min(Ntesty,[],2);    % true optimal design
      
        for temi = 1:nd    % design layer
            Naccurval = fun(desig(temi,:),Ncontext);
            NsampleCol = Naccurval+sd*randn(Nmc,n0);  %replications generator
            Nsample = mean(NsampleCol,2); tsd = sqrt(mean(var(NsampleCol,0,2))/n0);
            NGM1 = fitrgp(Ncontext,Nsample,'KernelFunction','exponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Nmy1(:,temi) = predict(NGM1,Ntestx);    % prediction of design temi at testx
            Nysd1 = (Nmy1(:,temi) - Ntesty(:,temi)).^2;
            Nms1(temi) = mean(Nysd1);
            
            NGM2 = fitrgp(Ncontext,Nsample,'KernelFunction','squaredexponential','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Nmy2(:,temi) = predict(NGM2,Ntestx);    % prediction of design temi at testx
            Nysd2 = (Nmy2(:,temi) - Ntesty(:,temi)).^2;
            Nms2(temi) = mean(Nysd2);
            
            NGM3 = fitrgp(Ncontext,Nsample,'KernelFunction','matern32','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Nmy3(:,temi) = predict(NGM3,Ntestx);    % prediction of design temi at testx
            Nysd3 = (Nmy3(:,temi) - Ntesty(:,temi)).^2;
            Nms3(temi) = mean(Nysd3);
            
            NGM4 = fitrgp(Ncontext,Nsample,'KernelFunction','matern52','FitMethod', ...
                'exact','PredictMethod','exact',... 
                'Sigma',tsd, 'ConstantSigma',true);    % GP for design temi
            Nmy4(:,temi) = predict(NGM4,Ntestx);    % prediction of design temi at testx
            Nysd4 = (Nmy4(:,temi) - Ntesty(:,temi)).^2;
            Nms4(temi) = mean(Nysd4);

        end
        NMmse1(mi) = max(Nms1);
        NMmse2(mi) = max(Nms2);
        NMmse3(mi) = max(Nms3);
        NMmse4(mi) = max(Nms4);
        
        [~,Napargmind1] = min(Nmy1,[],2);
        Nindtemp1 = (Napargmind1-1)*testa + ((1:testa)'); 
        Napargm1 = Ntesty(Nindtemp1);
        Napres1 = abs(Napargm1 - Nargm) > Del;    %whether false selection happen at testx
        NMpfs1(mi) = mean(Napres1);
        
        [~,Napargmind2] = min(Nmy2,[],2);
        Nindtemp2 = (Napargmind2-1)*testa + ((1:testa)'); 
        Napargm2 = Ntesty(Nindtemp2);
        Napres2 = abs(Napargm2 - Nargm) > Del;    %whether false selection happen at testx
        NMpfs2(mi) = mean(Napres2);
        
        [~,Napargmind3] = min(Nmy3,[],2);
        Nindtemp3 = (Napargmind3-1)*testa + ((1:testa)'); 
        Napargm3 = Ntesty(Nindtemp3);
        Napres3 = abs(Napargm3 - Nargm) > Del;    %whether false selection happen at testx
        NMpfs3(mi) = mean(Napres3);
        
        [~,Napargmind4] = min(Nmy4,[],2);
        Nindtemp4 = (Napargmind4-1)*testa + ((1:testa)'); 
        Napargm4 = Ntesty(Nindtemp4);
        Napres4 = abs(Napargm4 - Nargm) > Del;    %whether false selection happen at testx
        NMpfs4(mi) = mean(Napres4);
    end
    NPmse1(j) = mean(NMmse1); NPpfs1(j) = mean(NMpfs1);    %average of random covariates level
    NPmse2(j) = mean(NMmse2); NPpfs2(j) = mean(NMpfs2);    %average of random covariates level
    NPmse3(j) = mean(NMmse3); NPpfs3(j) = mean(NMpfs3);    %average of random covariates level
    NPmse4(j) = mean(NMmse4); NPpfs4(j) = mean(NMpfs4);    %average of random covariates level
end

end

