function [ UPmse1_all,UPmse2_all,UPmse3_all,...
    UPmse4_all,UPmse1,UPmse2,UPmse3,UPmse4,...
    pmU_mse_10,m_thres,UPnum1,UPnum2,UPnum3,UPnum4 ] ...
    = UniSQ_reg_self( m,desig,span1,span2,testa,Del,n0,cu,U,...
    thres1,thres )
%UniS returns the imse and ipfs under the uniform sampling
varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

rng(1000)
nm = length(m); nd = size(desig,1); xd = size(desig,2); sercost = cu*desig;
UPmse1_all = zeros(nm,1);     % number of covariates layer
UPmse2_all = zeros(nm,1);     % number of covariates layer
UPmse3_all = zeros(nm,1);     % number of covariates layer
UPmse4_all = zeros(nm,1);     % number of covariates layer

Ucontext_all = rand(max(m),1)*span1+span2;    % random covarites generator
Utestx = rand(testa,1)*span1+span2;   % test covariate generation for PFS
% Utesty = BestDesig(desig,Utestx,cu,U);    % exact testfun value
% Uargm = min(Utesty,[],2);    % true optimal design
Usample_all = zeros(max(m),nd); Usamplevar_all = zeros(max(m),nd);
for temi = 1:nd
    for si = 1:max(m)
        avrtime = QueueGG1(desig(temi),Ucontext_all(si),n0);
        Usample_all( si,temi ) = mean(min(avrtime + sercost(temi), U));
        Usamplevar_all( si,temi ) = var(min(avrtime + sercost(temi), U));
    end
end

num_subsets = 20*ones(length(m),1); num_subsets(end) = 1;
for j = 1:nm    % number of covariates layer
    j
    Umc = m(j);    
    UPmse1 = zeros(num_subsets(j),1);     % number of covariates layer
    UPmse2 = zeros(num_subsets(j),1);     % number of covariates layer
    UPmse3 = zeros(num_subsets(j),1);     % number of covariates layer
    UPmse4 = zeros(num_subsets(j),1);     % number of covariates layer
    for subseti = 1:num_subsets(j)
        UGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
        UGM2 = cell(nd,1); KernParac2 = cell(nd,1); coinvc2 = cell(nd,1);
        UGM3 = cell(nd,1); KernParac3 = cell(nd,1); coinvc3 = cell(nd,1);
        UGM4 = cell(nd,1); KernParac4 = cell(nd,1); coinvc4 = cell(nd,1);
        fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);
        fbasis2 = cell(nd,1); FSc2 = cell(nd,1); FSFc2 = cell(nd,1);
        fbasis3 = cell(nd,1); FSc3 = cell(nd,1); FSFc3 = cell(nd,1);
        fbasis4 = cell(nd,1); FSc4 = cell(nd,1); FSFc4 = cell(nd,1);
        Basf = ones(Umc,1);
        sample_ind = randperm(max(m), Umc);

        for temi = 1:nd    % design layer
            Ucontext = Ucontext_all(sample_ind);
            Usample = Usample_all(sample_ind,temi); 
            Usamplevar = Usamplevar_all(sample_ind,temi);
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
        kernind = 1;
        for temi = 1:nd
            KernParac1{temi} = [UGM1{temi}.theta,sqrt(UGM1{temi}.tausquared)];
            coinvc1{temi} = UGM1{temi}.Sigmainv;
            estbeta = UGM1{temi}.beta;
            fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis1{temi}(Ucontext))';
            FSc1{temi} = Fmat'*coinvc1{temi};
            FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
        end
        UPmse1(subseti) = mean(MSEEva_pred( Utestx,Ucontext,KernParac1,coinvc1,...
            FSFc1,FSc1,fbasis1,kernind,nd));

        kernind = 2;
        for temi = 1:nd
            KernParac2{temi} = [UGM2{temi}.theta,sqrt(UGM2{temi}.tausquared)];
            coinvc2{temi} = UGM2{temi}.Sigmainv;
            estbeta = UGM2{temi}.beta;
            fbasis2{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis2{temi}(Ucontext))';
            FSc2{temi} = Fmat'*coinvc2{temi};
            FSFc2{temi} = (FSc2{temi}*Fmat)^(-1);
        end
        UPmse2(subseti) = mean(MSEEva_pred(Utestx,Ucontext,KernParac2,coinvc2,...
            FSFc2,FSc2,fbasis2,kernind,nd));

        kernind = 3;
        for temi = 1:nd
            KernParac3{temi} = [UGM3{temi}.theta,sqrt(UGM3{temi}.tausquared)];
            coinvc3{temi} = UGM3{temi}.Sigmainv;
            estbeta = UGM3{temi}.beta;
            fbasis3{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis3{temi}(Ucontext))';
            FSc3{temi} = Fmat'*coinvc3{temi};
            FSFc3{temi} = (FSc3{temi}*Fmat)^(-1);
        end
        UPmse3(subseti) =  mean(MSEEva_pred(Utestx,Ucontext,KernParac3,coinvc3,...
            FSFc3,FSc3,fbasis3,kernind,nd));

        kernind = 4;
        for temi = 1:nd
            KernParac4{temi} = [UGM4{temi}.theta,sqrt(UGM4{temi}.tausquared)];
            coinvc4{temi} = UGM4{temi}.Sigmainv;
            estbeta = UGM4{temi}.beta;
            fbasis4{temi} = @(x) estbeta*ones(1,size(x,1));
            Fmat = (fbasis4{temi}(Ucontext))';
            FSc4{temi} = Fmat'*coinvc4{temi};
            FSFc4{temi} = (FSc4{temi}*Fmat)^(-1);
        end
        UPmse4(subseti) =  mean(MSEEva_pred(Utestx,Ucontext,KernParac4,coinvc4,...
            FSFc4,FSc4,fbasis4,kernind,nd));
    end
    UPmse1_all(j) = mean(UPmse1);     % number of covariates layer
    UPmse2_all(j) = mean(UPmse2);     % number of covariates layer
    UPmse3_all(j) = mean(UPmse3);     % number of covariates layer
    UPmse4_all(j) = mean(UPmse4);     % number of covariates layer

end

thres_logy1 = log(thres1);
thres_logy = log(thres);
pmU_mse_10 = zeros(4,2); % coefficients
px = log( m' );
m_thres = zeros(4,1);
% exponential:
t = 1;
py = log(UPmse1_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy1 - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% squared exponential:
t = 2;
py = log(UPmse2_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% matern 3/2:
t = 3;
py = log(UPmse3_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));
% matern 5/2:
t = 4;
py = log(UPmse4_all);
pmU_mse_10(t,:) = polyfit(px,py,1);
m_thres(t) = ceil(exp( (thres_logy - pmU_mse_10(t,2))/pmU_mse_10(t,1) ));

extram = m_thres - max(m);
Lm = 40;
kernind = 1;
if m_thres(kernind) > max(m)
[ UPmse1,UPnum1 ] = UniSQ_inner( desig,span1,span2,Lm,testa,Del,n0,cu,U,...
    Usample_all, Usamplevar_all, Ucontext_all, kernind, m_thres(kernind),thres1 );
else
    UPmse1 = UPmse1_all(end); UPnum1=0;
end
kernind = 2;
if m_thres(kernind) > max(m)
[ UPmse2,UPnum2 ] = UniSQ_inner( desig,span1,span2,Lm,testa,Del,n0,cu,U,...
    Usample_all, Usamplevar_all, Ucontext_all, kernind, m_thres(kernind),thres);
else
    UPmse2 = UPmse2_all(end); UPnum2=0;
end
kernind = 3;
if m_thres(kernind) > max(m)
[ UPmse3,UPnum3 ] = UniSQ_inner( desig,span1,span2,Lm,testa,Del,n0,cu,U,...
    Usample_all, Usamplevar_all, Ucontext_all, kernind, m_thres(kernind),thres);
else
    UPmse3 = UPmse3_all(end); UPnum3=0;
end
kernind = 4;
if m_thres(kernind) > max(m)
[ UPmse4,UPnum4 ] = UniSQ_inner( desig,span1,span2,Lm,testa,Del,n0,cu,U,...
    Usample_all, Usamplevar_all, Ucontext_all, kernind, m_thres(kernind),thres);
else
    UPmse4 = UPmse4_all(end); UPnum4=0;
end



end


function [ UPmse1,UPnum ] ...
    = UniSQ_inner( desig,span1,span2,Lm,testa,Del,n0,cu,U,...
    Usample1, Usamplevar1, Ucontext1, kernind, m_thres, vthres)

varliml = [0.001, 0.0007, 0.0006, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004,...
    0.0004, 0.0004];
varlimu = [0.15, 0.12, 0.09, 0.06, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02];

snum = length(Usample1); nd = size(desig,1); sercost = cu*desig;

UMmse1 = zeros(Lm,1); 
for mi = 1:Lm    % number of random covariates layer
    mi
    Utestx = rand(testa,1)*span1+span2;   % test covariate generation for PFS
    Utesty = BestDesig(desig,Utestx,cu,U);    % exact testfun value
    Uargm = min(Utesty,[],2);    % true optimal design
    Ums1 = zeros(nd,1); Umy1 = zeros(testa,nd);  
    Basf = ones(m_thres,1);
    Bpred = ones(testa,1);
    Ucontext = zeros(m_thres,1);   % random covarites generator
    Ucontext(1:snum) = Ucontext1;
    Ucontext((snum+1):end) = rand(m_thres-snum,1)*span1+span2; 
    
    UGM1 = cell(nd,1); KernParac1 = cell(nd,1); coinvc1 = cell(nd,1);
    fbasis1 = cell(nd,1); FSc1 = cell(nd,1); FSFc1 = cell(nd,1);

    for temi = 1:nd    % design layer
        Usample = zeros(m_thres,1);
        Usample(1:snum) = Usample1(:,temi);
        Usamplevar = zeros(m_thres,1);
        Usamplevar(1:snum) = Usamplevar1(:,temi);

        for si = (snum+1):m_thres
            avrtime = QueueGG1(desig(temi),Ucontext(si),n0);
            Usample( si ) = mean(min(avrtime + sercost(temi), U));
            Usamplevar( si ) = var(min(avrtime + sercost(temi), U));
        end
        estvar_poly = polyfit(Ucontext,Usamplevar,3);
        Usamplevar = max(min(polyval(estvar_poly,Ucontext), varlimu(temi)), varliml(temi))/n0;
        UGM1{temi} = SKfitnew(Ucontext, Usample, Basf, Usamplevar, kernind);
    
    end
    
    for temi = 1:nd
        KernParac1{temi} = [UGM1{temi}.theta,sqrt(UGM1{temi}.tausquared)];
        coinvc1{temi} = UGM1{temi}.Sigmainv;
        estbeta = UGM1{temi}.beta;
        fbasis1{temi} = @(x) estbeta*ones(1,size(x,1));
        Fmat = (fbasis1{temi}(Ucontext))';
        FSc1{temi} = Fmat'*coinvc1{temi};
        FSFc1{temi} = (FSc1{temi}*Fmat)^(-1);
    end
    UMmse1(mi) = mean(MSEEva_pred( Utestx,Ucontext,KernParac1,coinvc1,...
        FSFc1,FSc1,fbasis1,kernind,nd));
    

end

UPmse1 = mean(UMmse1);    % number of covariates layer
UPnum = UMmse1;




end










