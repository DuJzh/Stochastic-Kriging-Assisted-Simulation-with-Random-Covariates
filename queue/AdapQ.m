
clearvars,clc;

% MM1 Queue
fun = @QueueGG1; 
n0 = 10; n2 = 5; 
span1 = 4; span2 = 0.5;
desig = 6 + 0.3*(1:10); desig = desig'; 
nd = length(desig);  %location of design point
Lm = 100;    % number of replications for random covariates point

cu = 0.1; U = 2.5; 

m = 5*2.^(0:7); nm = length(m);    % number of randomly selected covariate
% hypothesis testing
intervnum = 10; tlim = max(m); tnum = 1:tlim; % from the start point to the end point
intdiv = (0:intervnum)/intervnum * span1 + span2;
crit = 16.92;
%%% PFS
testa = 1000;    % number of test covariate point for PFS
Del = 0.01;
%%%%%%% 
d = 1;  % context dim.
%%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
[ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniSQf( m,desig,span1,span2,Lm,testa,Del,n0,cu,U );
[ U5Pmse1,U5Ppfs1,U5Pmse2,U5Ppfs2,U5Pmse3,U5Ppfs3,U5Pmse4,U5Ppfs4 ] ...
    = UniSQf( m,desig,span1,span2,Lm,testa,Del,n2,cu,U );

%%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Truva = 3;
[ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnSQf( m,desig,span1,span2,Lm,testa,Del,n0,Truva,cu,U );
[ T5Pmse1,T5Ppfs1,T5Pmse2,T5Ppfs2,T5Pmse3,T5Ppfs3,T5Pmse4,T5Ppfs4 ] ...
    = TrnSQf( m,desig,span1,span2,Lm,testa,Del,n2,Truva,cu,U );




save testq__self;

funname = 'M/M/1 Queue';

figure
subplot(1,2,1)
hold on;
plot(m,U5Pmse2,'o-.')
plot(m,U5Pmse4,'d-.')
plot(m,U5Pmse3,'+-.')
plot(m,U5Pmse1,'*-.')
plot(m,UPmse2,'o-')
plot(m,UPmse4,'d-')
plot(m,UPmse3,'+-')
plot(m,UPmse1,'*-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Sq-Exp, n=5','Matern52, n=5','Matern32, n=5','Exp, n=5',...
    'Sq-Exp, n=10','Matern52, n=10','Matern32, n=10','Exp, n=10');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,2,2)
hold on;
plot(m,T5Pmse2,'o-.')
plot(m,T5Pmse4,'d-.')
plot(m,T5Pmse3,'+-.')
plot(m,T5Pmse1,'*-.')
plot(m,TPmse2,'o-')
plot(m,TPmse4,'d-')
plot(m,TPmse3,'+-')
plot(m,TPmse1,'*-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Sq-Exp, n=5','Matern52, n=5','Matern32, n=5','Exp, n=5',...
    'Sq-Exp, n=10','Matern52, n=10','Matern32, n=10','Exp, n=10');
% set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')



figure
subplot(1,2,1)
hold on;
plot(m,U5Ppfs2,'o-.')
plot(m,U5Ppfs4,'d-.')
plot(m,U5Ppfs3,'+-.')
plot(m,U5Ppfs1,'*-.')
plot(m,UPpfs2,'o-')
plot(m,UPpfs4,'d-')
plot(m,UPpfs3,'+-')
plot(m,UPpfs1,'*-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Sq-Exp, n=5','Matern52, n=5','Matern32, n=5','Exp, n=5',...
    'Sq-Exp, n=10','Matern52, n=10','Matern32, n=10','Exp, n=10');
set(gca,'XTick', m, 'xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,2,2)
hold on;
plot(m,T5Ppfs2,'o-.')
plot(m,T5Ppfs4,'d-.')
plot(m,T5Ppfs3,'+-.')
plot(m,T5Ppfs1,'*-.')
plot(m,TPpfs2,'o-')
plot(m,TPpfs4,'d-')
plot(m,TPpfs3,'+-')
plot(m,TPpfs1,'*-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Sq-Exp, n=5','Matern52, n=5','Matern32, n=5','Exp, n=5',...
    'Sq-Exp, n=10','Matern52, n=10','Matern32, n=10','Exp, n=10');
% set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')






