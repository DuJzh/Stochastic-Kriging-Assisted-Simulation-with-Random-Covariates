function [ TPmse1_all,TPmse2_all,TPmse3_all,...
    TPmse4_all,TPmse1,TPmse2,TPmse3,TPmse4,...
    pmU_mse_10,m_thres,TPnum1,TPnum2,TPnum3,TPnum4 ] ...
    = TrnSQ_reg_self( m,desig,span1,span2,testa,Del,n0,Truva,cu,U,...
    thres1,thres )
%TrnS returns the imse and ipfs under the truncated normal sampling
varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

rng(1000)
nm = length(m); nd = size(desig,1); xd = size(desig,2); sercost = cu*desig;
Tmu = span2+span1/2;
pd = makedist('Normal','mu',Tmu,'sigma',Truva); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);

TPmse1_all = zeros(nm,1);  % number of covariates layer
TPmse2_all = zeros(nm,1);  % number of covariates layer
TPmse3_all = zeros(nm,1);  % number of covariates layer
TPmse4_all = zeros(nm,1);  % number of covariates layer

Tcontext_all = random(Tpd,max(m),1);     % random covarites generator
Ttestx = random(Tpd,testa,1);            % test covariate generation for PFS
% Ttesty = BestDesig(desig,Ttestx,cu,U);    % exact testfun value
% Targm = min(Ttesty,[],2);    % true optimal design
Tsample_all = zeros(max(m),nd); Tsamplevar_all = zeros(max(m),nd);
for temi = 1:nd
    for si = 1:max(m)
        avrtime = QueueGG1(desig(temi),Tcontext_all(si),n0);
        Tsample_all( si,temi ) = mean(min(avrtime + sercost(temi), U));
        Tsamplevar_all( si,temi ) = var(min(avrtime + sercost(temi), U));
    end
end

num_subsets = 20*ones(length(m),1); num_subsets(end) = 1;
for j = 1:nm    % number of covariates layer
    j
    Tmc = m(j);    
    TPmse1 = zeros(num_subsets(j),1);     % number of covariates layer
    TPmse2 = zeros(num_subsets(j),1);     % number of covariates layer
    TPmse3 = zeros(num_subsets(j),1);     % number of covariates layer
    TPmse4 = zeros(num_subsets(j),1);     % number of covariates layer
    for subseti = 1:num_subsets(j)
        TGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
        TGM2 = cell(nd,1); KernParac2 = cell(nd,1); coinvc2 = cell(nd,1);
        TGM3 = cell(nd,1); KernParac3 = cell(nd,1); coinvc3 = cell(nd,1);
        TGM4 = cell(nd,1); KernParac4 = cell(nd,1); coinvc4 = cell(nd,1);
        fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);
        fbasis2 = cell(nd,1); FSc2 = cell(nd,1); FSFc2 = cell(nd,1);
        fbasis3 = cell(nd,1); FSc3 = cell(nd,1); FSFc3 = cell(nd,1);
        fbasis4 = cell(nd,1); FSc4 = cell(nd,1); FSFc4 = cell(nd,1);
        Basf = ones(Tmc,1);
        sample_ind = randperm(max(m), Tmc);

        for temi = 1:nd    % design layer
            Tcontext = Tcontext_all(sample_ind);
            Tsample = Tsample_all(sample_ind,temi); 
            Tsamplevar = Tsamplevar_all(sample_ind,temi);
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
        kernind = 1;
        for temi = 1:nd
            KernParac1{temi} = [TGM1{temi}.theta,sqrt(TGM1{temi}.tausquared)];
            coinvc1{temi} = TGM1{temi}.Sigmainv;
            estbeta = TGM1{temi}.beta;
            fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis1{temi}(Tcontext))';
            FSc1{temi} = Fmat'*coinvc1{temi};
            FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
        end
        TPmse1(subseti) = mean(MSEEva_pred( Ttestx,Tcontext,KernParac1,coinvc1,...
            FSFc1,FSc1,fbasis1,kernind,nd));

        kernind = 2;
        for temi = 1:nd
            KernParac2{temi} = [TGM2{temi}.theta,sqrt(TGM2{temi}.tausquared)];
            coinvc2{temi} = TGM2{temi}.Sigmainv;
            estbeta = TGM2{temi}.beta;
            fbasis2{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis2{temi}(Tcontext))';
            FSc2{temi} = Fmat'*coinvc2{temi};
            FSFc2{temi} = (FSc2{temi}*Fmat)^(-1);
        end
        TPmse2(subseti) = mean(MSEEva_pred(Ttestx,Tcontext,KernParac2,coinvc2,...
            FSFc2,FSc2,fbasis2,kernind,nd));

        kernind = 3;
        for temi = 1:nd
            KernParac3{temi} = [TGM3{temi}.theta,sqrt(TGM3{temi}.tausquared)];
            coinvc3{temi} = TGM3{temi}.Sigmainv;
            estbeta = TGM3{temi}.beta;
            fbasis3{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis3{temi}(Tcontext))';
            FSc3{temi} = Fmat'*coinvc3{temi};
            FSFc3{temi} = (FSc3{temi}*Fmat)^(-1);
        end
        TPmse3(subseti) =  mean(MSEEva_pred(Ttestx,Tcontext,KernParac3,coinvc3,...
            FSFc3,FSc3,fbasis3,kernind,nd));

        kernind = 4;
        for temi = 1:nd
            KernParac4{temi} = [TGM4{temi}.theta,sqrt(TGM4{temi}.tausquared)];
            coinvc4{temi} = TGM4{temi}.Sigmainv;
            estbeta = TGM4{temi}.beta;
            fbasis4{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis4{temi}(Tcontext))';
            FSc4{temi} = Fmat'*coinvc4{temi};
            FSFc4{temi} = (FSc4{temi}*Fmat)^(-1);
        end
        TPmse4(subseti) =  mean(MSEEva_pred(Ttestx,Tcontext,KernParac4,coinvc4,...
            FSFc4,FSc4,fbasis4,kernind,nd));
    end
    TPmse1_all(j) = mean(TPmse1); % number of covariates layer
    TPmse2_all(j) = mean(TPmse2); % number of covariates layer
    TPmse3_all(j) = mean(TPmse3); % number of covariates layer
    TPmse4_all(j) = mean(TPmse4); % number of covariates layer
    
end

thres_logy1 = log(thres1);
thres_logy = log(thres);
pmU_mse_10 = zeros(4,2); % coefficients
px = log( m' );
m_thres = zeros(4,1);
% exponential:
t = 1;
py = log(TPmse1_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy1 - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% squared exponential:
t = 2;
py = log(TPmse2_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% matern 3/2:
t = 3;
py = log(TPmse3_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% matern 5/2:
t = 4;
py = log(TPmse4_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));

extram = m_thres - max(m);
Lm = 40;
kernind = 1;
if m_thres(kernind) > max(m)
[ TPmse1,TPnum1 ] = TrnSQ_inner( desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U,...
    Tsample_all, Tsamplevar_all, Tcontext_all, kernind, m_thres(kernind),thres1);
else
    TPmse1 = TPmse1_all(end); TPnum1=0;
end
kernind = 2;
if m_thres(kernind) > max(m)
[ TPmse2,TPnum2 ] = TrnSQ_inner( desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U,...
    Tsample_all, Tsamplevar_all, Tcontext_all, kernind, m_thres(kernind),thres);
else
    TPmse2 = TPmse2_all(end); TPnum2=0;
end
kernind = 3;
if m_thres(kernind) > max(m)
[ TPmse3,TPnum3 ] = TrnSQ_inner( desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U,...
    Tsample_all, Tsamplevar_all, Tcontext_all, kernind, m_thres(kernind),thres);
else
    TPmse3 = TPmse3_all(end); TPnum3=0;
end
kernind = 4;
if m_thres(kernind) > max(m)
[ TPmse4,TPnum4 ] = TrnSQ_inner( desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U,...
    Tsample_all, Tsamplevar_all, Tcontext_all, kernind, m_thres(kernind),thres);
else
    TPmse4 = TPmse4_all(end); TPnum4=0;
end



end



function [ TPmse1,TPnum ] ...
    = TrnSQ_inner( desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U,...
    Tsample1, Tsamplevar1, Tcontext1, kernind, m_thres, vthres)

varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

snum = length(Tsample1); nd = size(desig,1); sercost = cu*desig;

Tmu = span2+span1/2;
pd = makedist('Normal','mu',Tmu,'sigma',Truva); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);

TMmse1 = zeros(Lm,1); 
for mi = 1:Lm    % number of random covariates layer
    mi
    Ttestx = random(Tpd,testa,1);      % test covariate generation for PFS
    Ttesty = BestDesig(desig,Ttestx,cu,U);    % exact testfun value
    Targm = min(Ttesty,[],2);    % true optimal design
    Tms1 = zeros(nd,1); Tmy1 = zeros(testa,nd);  
    Basf = ones(m_thres,1);
    Bpred = ones(testa,1);
    Tcontext = zeros(m_thres,1);   % random covarites generator
    Tcontext(1:snum) = Tcontext1;
    Tcontext((snum+1):end) = random(Tpd,(m_thres-snum),1); 
    
    TGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
    fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);

    for temi = 1:nd    % design layer
        Tsample = zeros(m_thres,1);
        Tsample(1:snum) = Tsample1(:,temi);
        Tsamplevar = zeros(m_thres,1);
        Tsamplevar(1:snum) = Tsamplevar1(:,temi);

        for si = (snum+1):m_thres
            avrtime = QueueGG1(desig(temi),Tcontext(si),n0);
            Tsample( si ) = mean(min(avrtime + sercost(temi), U));
            Tsamplevar( si ) = var(min(avrtime + sercost(temi), U));
        end
        estvar_poly = polyfit(Tcontext,Tsamplevar,3);
        Tsamplevar = max(min(polyval(estvar_poly,Tcontext), varlimu(temi)), varliml(temi))/n0;
        TGM1{temi} = SKfitnew(Tcontext, Tsample, Basf, Tsamplevar, kernind);
    
    end
    
    for temi = 1:nd
        KernParac1{temi} = [TGM1{temi}.theta,sqrt(TGM1{temi}.tausquared)];
        coinvc1{temi} = TGM1{temi}.Sigmainv;
        estbeta = TGM1{temi}.beta;
        fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
        Fmat = (fbasis1{temi}(Tcontext))';
        FSc1{temi} = Fmat'*coinvc1{temi};
        FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
    end
    TMmse1(mi) = mean(MSEEva_pred( Ttestx,Tcontext,KernParac1,coinvc1,...
        FSFc1,FSc1,fbasis1,kernind,nd));
end

TPmse1 = mean(TMmse1);   % number of covariates layer
TPnum = TMmse1;




end









