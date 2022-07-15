
clearvars,clc;

% MM1 Queue
fun = @QueueGG1; 
n0 = 10; n2 = 5; 
span1 = 4; span2 = 0.5;
desig = 6 + 0.3*(1:10); desig = desig'; 
nd = length(desig);  %location of design point
Lm = 100;    % number of replications for random covariates point

cu = 0.1; U = 2.5; 

m = [10,15,23,35,53,80]; nm = length(m);    % number of randomly selected covariate
% hypothesis testing
intervnum = 10; tlim = max(m); tnum = 1:tlim; % from the start point to the end point
intdiv = (0:intervnum)/intervnum * span1 + span2;
crit = 16.92;
%%% PFS
testa = 40000;    % number of test covariate point for PFS
Del = 0.01;
%%%%%%%
% lose the error tolerance in finding the optimal point
options = optimoptions('fmincon','OptimalityTolerance',0.005,...
    'StepTolerance',0.005,'Display','off');
optseed = 10;  
d = 1;  % context dim.
% %%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
thres1 = 0.0002;  % 2
thres2 = 0.000075;  %075
[ UPmse1_all,UPmse2_all,UPmse3_all,...
    UPmse4_all,pUPmse1,pUPmse2,pUPmse3,pUPmse4,pmU_mse,m_thresU,...
    UPnum1,UPnum2,UPnum3,UPnum4] ...
    = UniSQ_reg_self( m,desig,span1,span2,testa,Del,n0,cu,U,thres1,thres2 );

% %%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Truva = 3;
[ TPmse1_all,TPmse2_all,TPmse3_all,...
    TPmse4_all,pTPmse1,pTPmse2,pTPmse3,pTPmse4,pmT_mse,m_thresT,...
    TPnum1,TPnum2,TPnum3,TPnum4] ...
    = TrnSQ_reg_self( m,desig,span1,span2,testa,Del,n0,Truva,cu,U,thres1,thres2 );

save testq_reg;

funname = 'M/M/1 Queue';

%%%%%%%%%%%

figure
subplot(1,2,1)
plot(m,UPmse1,'*-')
hold on;
plot(m,UPmse2,'o-')
plot(m,UPmse3,'+-')
plot(m,UPmse4,'d-')
plot(m,U5Pmse1,'*-.')
plot(m,U5Pmse2,'o-.')
plot(m,U5Pmse3,'+-.')
plot(m,U5Pmse4,'d-.')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Exp, n=10','Sq-Exp, n=10','Matern32, n=10','Matern52, n=10',...
    'Exp, n=5','Sq-Exp, n=5','Matern32, n=5','Matern52, n=5');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,2,2)
plot(m,TPmse1,'*-')
hold on;
plot(m,TPmse2,'o-')
plot(m,TPmse3,'+-')
plot(m,TPmse4,'d-')
plot(m,T5Pmse1,'*-.')
plot(m,T5Pmse2,'o-.')
plot(m,T5Pmse3,'+-.')
plot(m,T5Pmse4,'d-.')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Exp, n=10','Sq-Exp, n=10','Matern32, n=10','Matern52, n=10',...
    'Exp, n=5','Sq-Exp, n=5','Matern32, n=5','Matern52, n=5');
% set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')



figure
subplot(1,2,1)
plot(m,UPpfs1,'*-')
hold on;
plot(m,UPpfs2,'o-')
plot(m,UPpfs3,'+-')
plot(m,UPpfs4,'d-')
plot(m,U5Ppfs1,'*-.')
plot(m,U5Ppfs2,'o-.')
plot(m,U5Ppfs3,'+-.')
plot(m,U5Ppfs4,'d-.')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Exp, n=10','Sq-Exp, n=10','Matern32, n=10','Matern52, n=10',...
    'Exp, n=5','Sq-Exp, n=5','Matern32, n=5','Matern52, n=5');
set(gca,'XTick', m, 'xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,2,2)
plot(m,TPpfs1,'*-')
hold on;
plot(m,TPpfs2,'o-')
plot(m,TPpfs3,'+-')
plot(m,TPpfs4,'d-')
plot(m,T5Ppfs1,'*-.')
plot(m,T5Ppfs2,'o-.')
plot(m,T5Ppfs3,'+-.')
plot(m,T5Ppfs4,'d-.')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Exp, n=10','Sq-Exp, n=10','Matern32, n=10','Matern52, n=10',...
    'Exp, n=5','Sq-Exp, n=5','Matern32, n=5','Matern52, n=5');
% set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')









