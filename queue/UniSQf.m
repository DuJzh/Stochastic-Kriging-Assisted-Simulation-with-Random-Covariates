function [ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniSQf( m,desig,span1,span2,Lm,testa,Del,n0,cu,U )
%UniS returns the imse and ipfs under the uniform sampling

% lower and upper boound of the variance for each design
varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

nm = length(m); nd = size(desig,1); xd = size(desig,2); sercost = cu*desig;
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
        ts = sercost; tdesig = desig;
        Ucontext = rand(Umc,1)*span1+span2;    % random covarites generator  
        Utestx = rand(testa,1)*span1+span2;   % test covariate generation for PFS
        Basf = ones(Umc,1);
        Bpred = ones(testa,1);

        UGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
        UGM2 = cell(nd,1); KernParac2 = cell(nd,1); coinvc2 = cell(nd,1);
        UGM3 = cell(nd,1); KernParac3 = cell(nd,1); coinvc3 = cell(nd,1);
        UGM4 = cell(nd,1); KernParac4 = cell(nd,1); coinvc4 = cell(nd,1);
        fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);
        fbasis2 = cell(nd,1); FSc2 = cell(nd,1); FSFc2 = cell(nd,1);
        fbasis3 = cell(nd,1); FSc3 = cell(nd,1); FSFc3 = cell(nd,1);
        fbasis4 = cell(nd,1); FSc4 = cell(nd,1); FSFc4 = cell(nd,1);

        for temi = 1:nd    % design layer
            Usample = zeros(Umc,1); Usamplevar = zeros(Umc,1);
            for si = 1:Umc
                avrtime = QueueGG1(desig(temi),Ucontext(si),n0);
                Usample( si ) = mean(min(avrtime + sercost(temi), U));
                Usamplevar( si ) = var(min(avrtime + sercost(temi), U));
            end
            estvar_poly = polyfit(Ucontext,Usamplevar,3);
            Usamplevar = max(min(polyval(estvar_poly,Ucontext), varlimu(temi)), varliml(temi))/n0;
            
            kernind = 1;
            UGM1{temi} = SKfitnew(Ucontext, Usample, Basf, Usamplevar, kernind);

            kernind = 2;
            UGM2{temi} = SKfitnew(Ucontext, Usample, Basf, Usamplevar, kernind);

            kernind = 3;
            UGM3{temi} =  SKfitnew(Ucontext, Usample, Basf, Usamplevar, kernind);

            kernind = 4;
            UGM4{temi} = SKfitnew(Ucontext, Usample, Basf, Usamplevar, kernind);
        end
        fi1 = cell(nd,1); kernind = 1;
        for temi = 1:nd
            KernParac1{temi} = [UGM1{temi}.theta,sqrt(UGM1{temi}.tausquared)];
            coinvc1{temi} = UGM1{temi}.Sigmainv;
            estbeta = UGM1{temi}.beta;
            fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis1{temi}(Ucontext))';
            FSc1{temi} = Fmat'*coinvc1{temi};
            FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
            fi1{temi} = @(testx) MSEEva( testx,Ucontext,KernParac1{temi},...
                coinvc1{temi},FSFc1{temi},FSc1{temi},fbasis1{temi},kernind );
        end
        UMpfs1(mi) = mean(PFSEva3( Utestx,UGM1,nd,Del,fi1,xd,Bpred ));
        UMmse1(mi) = mean(MSEEva_pred( Utestx,Ucontext,KernParac1,coinvc1,...
            FSFc1,FSc1,fbasis1,kernind,nd));

        fi2 = cell(nd,1); kernind = 2;
        for temi = 1:nd
            KernParac2{temi} = [UGM2{temi}.theta,sqrt(UGM2{temi}.tausquared)];
            coinvc2{temi} = UGM2{temi}.Sigmainv;
            estbeta = UGM2{temi}.beta;
            fbasis2{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis2{temi}(Ucontext))';
            FSc2{temi} = Fmat'*coinvc2{temi};
            FSFc2{temi} = (FSc2{temi}*Fmat)^(-1);
            fi2{temi} = @(testx) MSEEva( testx,Ucontext,KernParac2{temi},...
                coinvc2{temi},FSFc2{temi},FSc2{temi},fbasis2{temi},kernind );
        end
        UMpfs2(mi) = mean(PFSEva3( Utestx,UGM2,nd,Del,fi2,xd,Bpred ));
        UMmse2(mi) = mean(MSEEva_pred(Utestx,Ucontext,KernParac2,coinvc2,...
            FSFc2,FSc2,fbasis2,kernind,nd));

        fi3 = cell(nd,1); kernind = 3;
        for temi = 1:nd
            KernParac3{temi} = [UGM3{temi}.theta,sqrt(UGM3{temi}.tausquared)];
            coinvc3{temi} = UGM3{temi}.Sigmainv;
            estbeta = UGM3{temi}.beta;
            fbasis3{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis3{temi}(Ucontext))';
            FSc3{temi} = Fmat'*coinvc3{temi};
            FSFc3{temi} = (FSc3{temi}*Fmat)^(-1);
            fi3{temi} = @(testx) MSEEva( testx,Ucontext,KernParac3{temi},...
                coinvc3{temi},FSFc3{temi},FSc3{temi},fbasis3{temi},kernind );
        end
        UMpfs3(mi) = mean(PFSEva3( Utestx,UGM3,nd,Del,fi3,xd,Bpred ));
        UMmse3(mi) =  mean(MSEEva_pred(Utestx,Ucontext,KernParac3,coinvc3,...
            FSFc3,FSc3,fbasis3,kernind,nd));

        fi4 = cell(nd,1); kernind = 4;
        for temi = 1:nd
            KernParac4{temi} = [UGM4{temi}.theta,sqrt(UGM4{temi}.tausquared)];
            coinvc4{temi} = UGM4{temi}.Sigmainv;
            estbeta = UGM4{temi}.beta;
            fbasis4{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis4{temi}(Ucontext))';
            FSc4{temi} = Fmat'*coinvc4{temi};
            FSFc4{temi} = (FSc4{temi}*Fmat)^(-1);
            fi4{temi} = @(testx) MSEEva( testx,Ucontext,KernParac4{temi},...
                coinvc4{temi},FSFc4{temi},FSc4{temi},fbasis4{temi},kernind );
        end
        UMpfs4(mi) = mean(PFSEva3( Utestx,UGM4,nd,Del,fi4,xd,Bpred ));
        UMmse4(mi) =  mean(MSEEva_pred(Utestx,Ucontext,KernParac4,coinvc4,...
            FSFc4,FSc4,fbasis4,kernind,nd));
    end
    UPmse1(j) = mean(UMmse1); UPpfs1(j) = mean(UMpfs1);    %average of random covariates level
    UPmse2(j) = mean(UMmse2); UPpfs2(j) = mean(UMpfs2);    %average of random covariates level
    UPmse3(j) = mean(UMmse3); UPpfs3(j) = mean(UMpfs3);    %average of random covariates level
    UPmse4(j) = mean(UMmse4); UPpfs4(j) = mean(UMpfs4);    %average of random covariates level
end

end



