function [ EPmse1,EPpfs1,Econshortrec,TT ] ...
    = AdaS2D( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d)
%AdaS returns the imse and ipfs under the adaptive sampling

if kernind == 1  %exponential
    kername = 'exponential';
elseif kernind == 2  % squared-exponential
    kername = 'squaredexponential';
elseif kernind == 3  % matern 3/2
    kername = 'matern32';
elseif kernind == 4  % matern 5/2
    kername = 'matern52';
end


nm = length(m); nd = length(desig);
EMmse1 = zeros(nm,Lm); EMpfs1 = zeros(nm,Lm);    % test covariate point layer
interv = [span2, span2+span1]; Econshortrec = []; TT = [];
lb = interv(1)*ones(1,d); ub = interv(2)*ones(1,d);
for mi = 1:Lm    % number of random covariates layer
    mi
    tdesig = desig; tEtesty = Etesty; tm = m; tfun = fun; TT1 = [];
    EGM = cell(nd,1); KernParac = cell(nd,1); coinvc = cell(nd,1); 
    fbasis = cell(nd,1); FSc = cell(nd,1); FSFc = cell(nd,1);
    tvan = zeros(nd,1);
    Econshort = (1/2*span1+span2)*ones(1,d);
    Esample = cell(nd,1);
    for temi = 1:nd
        Eaccurval = fun(desig(temi,:),Econshort);
        Esampletrain = Eaccurval + sd*randn(n0,1);    %replications generator
        Esample{temi} = [Esample{temi}; mean(Esampletrain)];  % sample
        tvan(temi) = var(Esampletrain);  % var est
    end
    nextm = rand(1,d)*span1+span2; 
    Econshort = [Econshort;nextm]; dm = [m(1)-1,diff(m)]; snum = 1;
    for mm = 1:nm
        Ems1 = zeros(nd,1); Emy1 = zeros(testa,nd);  
        for Ei = 1:dm(mm)
            for temi = 1:nd    % design layer
                Eaccurval = fun(desig(temi,:),nextm);
                Esampletrain = Eaccurval + sd*randn(n0,1);    %replications generator
                Esample{temi} = [Esample{temi}; mean(Esampletrain)];  % sample
                tvan(temi) = (tvan(temi)*snum + var(Esampletrain))/(snum+1); % var est
                tsdn = sqrt(tvan(temi)/n0);
                EGM{temi} = fitrgp(Econshort,Esample{temi},'KernelFunction',...
                    kername,'FitMethod', ...
                    'exact','PredictMethod','exact',... 
                    'Sigma',tsdn, 'ConstantSigma',true);    % GP for design temi
            end
            snum = snum+1;
            tic
            for temi = 1:nd
                KernParac{temi} = EGM{temi}.KernelInformation.KernelParameters;
                estsd = EGM{temi}.Sigma;
                coinvc{temi} = (Kernval(Econshort,Econshort,KernParac{temi},kernind) ...
                    + estsd^2*eye(snum))^(-1);
                estbeta = EGM{temi}.Beta;
                fbasis{temi} = @(x) estbeta*ones(1,size(x,1));
                Fmat = (fbasis{temi}(Econshort))';
                FSc{temi} = Fmat'*coinvc{temi};
                FSFc{temi} = (FSc{temi}*Fmat)^(-1);
            end
            mmse = @(x) (-1)*MSEEva2D( x,Econshort,KernParac,coinvc,...
                FSFc,FSc,fbasis,kernind,nd );
            nextmm = []; tempy = []; 
            slist = rand(optseed,d)*span1+span2;  % different initial points
            for ni = 1:length(slist)  % nextmm contains the local minimum
                [tempmm,tempyy] = fmincon(mmse,slist(ni,:),[],[],[],[],...
                    lb,ub,[],options);
                nextmm = [nextmm;tempmm];
                tempy = [tempy;tempyy];
            end
            [~,indd] = min(tempy);  % the best local minimum
            TT1 = [TT1;toc];
            nextm = nextmm(indd,:); 
            Econshort = [Econshort;nextm];
        end
        for temi = 1:nd    % design layer
            Emy1(:,temi) = predict(EGM{temi},Etestx);    % prediction of design temi at testx
            Eysd1 = (Emy1(:,temi) - Etesty(:,temi)).^2;
            Ems1(temi) = mean(Eysd1);
        end
        EMmse1(mm,mi) = max(Ems1);
        [~,Eapargmind1] = min(Emy1,[],2);
        Eindtemp1 = (Eapargmind1-1)*testa + ((1:testa)'); 
        Eapargm1 = Etesty(Eindtemp1);
        Eapres1 = abs(Eapargm1 - Eargm) > Del;    %whether false selection happen at testx
        EMpfs1(mm,mi) = mean(Eapres1);
    end
    Econshortrec = [Econshortrec,Econshort];
    TT = [TT,TT1];
end
EPmse1 = mean(EMmse1,2); EPpfs1 = mean(EMpfs1,2);    %average of random covariates level

end

