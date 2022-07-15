function [ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnSQf( m,desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U )
%TrnS returns the imse and ipfs under the truncated normal sampling

% lower and upper boound of the variance for each design
varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

nm = length(m); nd = size(desig,1); xd = size(desig,2); sercost = cu*desig;
Tmu = span2+span1/2;
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
        ts = sercost; tdesig = desig;
        Tcontext = random(Tpd,Tmc,1);      % random covarites generator
        Ttestx = random(Tpd,testa,1);         % test covariate generation for PFS
        Basf = ones(Tmc,1);
        Bpred = ones(testa,1);

        TGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
        TGM2 = cell(nd,1); KernParac2 = cell(nd,1); coinvc2 = cell(nd,1);
        TGM3 = cell(nd,1); KernParac3 = cell(nd,1); coinvc3 = cell(nd,1);
        TGM4 = cell(nd,1); KernParac4 = cell(nd,1); coinvc4 = cell(nd,1);
        fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);
        fbasis2 = cell(nd,1); FSc2 = cell(nd,1); FSFc2 = cell(nd,1);
        fbasis3 = cell(nd,1); FSc3 = cell(nd,1); FSFc3 = cell(nd,1);
        fbasis4 = cell(nd,1); FSc4 = cell(nd,1); FSFc4 = cell(nd,1);

        for temi = 1:nd    % design layer
            Tsample = zeros(Tmc,1); Tsamplevar = zeros(Tmc,1);
            for si = 1:Tmc
                avrtime = QueueGG1(desig(temi),Tcontext(si),n0);
                Tsample( si ) = mean(min(avrtime + sercost(temi), U));
                Tsamplevar( si ) = var(min(avrtime + sercost(temi), U));
            end  
            estvar_poly = polyfit(Tcontext,Tsamplevar,3);
            Tsamplevar = max(min(polyval(estvar_poly,Tcontext), varlimu(temi)), varliml(temi))/n0;
            kernind = 1;
            TGM1{temi} = SKfitnew(Tcontext, Tsample, Basf, Tsamplevar, kernind);

            kernind = 2;
            TGM2{temi} = SKfitnew(Tcontext, Tsample, Basf, Tsamplevar, kernind);

            kernind = 3;
            TGM3{temi} = SKfitnew(Tcontext, Tsample, Basf, Tsamplevar, kernind);

            kernind = 4;
            TGM4{temi} = SKfitnew(Tcontext, Tsample, Basf, Tsamplevar, kernind);
        end
        fi1 = cell(nd,1);  kernind = 1;
        for temi = 1:nd
            KernParac1{temi} = [TGM1{temi}.theta,sqrt(TGM1{temi}.tausquared)];
            coinvc1{temi} = TGM1{temi}.Sigmainv;
            estbeta = TGM1{temi}.beta;
            fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis1{temi}(Tcontext))';
            FSc1{temi} = Fmat'*coinvc1{temi};
            FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
            fi1{temi} = @(testx) MSEEva( testx,Tcontext,KernParac1{temi},...
                coinvc1{temi},FSFc1{temi},FSc1{temi},fbasis1{temi},kernind );
        end
        TMpfs1(mi) = mean(PFSEva3( Ttestx,TGM1,nd,Del,fi1,xd,Bpred ));
        TMmse1(mi) = mean(MSEEva_pred( Ttestx,Tcontext,KernParac1,coinvc1,...
            FSFc1,FSc1,fbasis1,kernind,nd));

        fi2 = cell(nd,1); kernind = 2;
        for temi = 1:nd
            KernParac2{temi} = [TGM2{temi}.theta,sqrt(TGM2{temi}.tausquared)];
            coinvc2{temi} = TGM2{temi}.Sigmainv;
            estbeta = TGM2{temi}.beta;
            fbasis2{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis2{temi}(Tcontext))';
            FSc2{temi} = Fmat'*coinvc2{temi};
            FSFc2{temi} = (FSc2{temi}*Fmat)^(-1);
            fi2{temi} = @(testx) MSEEva( testx,Tcontext,KernParac2{temi},...
                coinvc2{temi},FSFc2{temi},FSc2{temi},fbasis2{temi},kernind );
        end
        TMpfs2(mi) = mean(PFSEva3( Ttestx,TGM2,nd,Del,fi2,xd,Bpred ));
        TMmse2(mi) = mean(MSEEva_pred(Ttestx,Tcontext,KernParac2,coinvc2,...
            FSFc2,FSc2,fbasis2,kernind,nd));

        fi3 = cell(nd,1); kernind = 3;
        for temi = 1:nd
            KernParac3{temi} = [TGM3{temi}.theta,sqrt(TGM3{temi}.tausquared)];
            coinvc3{temi} = TGM3{temi}.Sigmainv;
            estbeta = TGM3{temi}.beta;
            fbasis3{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis3{temi}(Tcontext))';
            FSc3{temi} = Fmat'*coinvc3{temi};
            FSFc3{temi} = (FSc3{temi}*Fmat)^(-1);
            fi3{temi} = @(testx) MSEEva( testx,Tcontext,KernParac3{temi},...
                coinvc3{temi},FSFc3{temi},FSc3{temi},fbasis3{temi},kernind );
        end
        TMpfs3(mi) = mean(PFSEva3( Ttestx,TGM3,nd,Del,fi3,xd,Bpred ));
        TMmse3(mi) =  mean(MSEEva_pred(Ttestx,Tcontext,KernParac3,coinvc3,...
            FSFc3,FSc3,fbasis3,kernind,nd));

        fi4 = cell(nd,1); kernind = 4;
        for temi = 1:nd
            KernParac4{temi} = [TGM4{temi}.theta,sqrt(TGM4{temi}.tausquared)];
            coinvc4{temi} = TGM4{temi}.Sigmainv;
            estbeta = TGM4{temi}.beta;
            fbasis4{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis4{temi}(Tcontext))';
            FSc4{temi} = Fmat'*coinvc4{temi};
            FSFc4{temi} = (FSc4{temi}*Fmat)^(-1);
            fi4{temi} = @(testx) MSEEva( testx,Tcontext,KernParac4{temi},...
                coinvc4{temi},FSFc4{temi},FSc4{temi},fbasis4{temi},kernind );
        end
        TMpfs4(mi) = mean(PFSEva3( Ttestx,TGM4,nd,Del,fi4,xd,Bpred ));
        TMmse4(mi) =  mean(MSEEva_pred(Ttestx,Tcontext,KernParac4,coinvc4,...
            FSFc4,FSc4,fbasis4,kernind,nd));
    
    end
    TPmse1(j) = mean(TMmse1); TPpfs1(j) = mean(TMpfs1);    %average of random covariates level
    TPmse2(j) = mean(TMmse2); TPpfs2(j) = mean(TMpfs2);    %average of random covariates level
    TPmse3(j) = mean(TMmse3); TPpfs3(j) = mean(TMpfs3);    %average of random covariates level
    TPmse4(j) = mean(TMmse4); TPpfs4(j) = mean(TMpfs4);    %average of random covariates level
end

end









