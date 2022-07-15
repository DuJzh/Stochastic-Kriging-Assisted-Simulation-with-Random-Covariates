
clearvars,clc;

% Griewangk's function: testfun2
fun = @testfun2;
desig = 1:10; desig = desig'; d = 10;  % context dim.  desig=1:4;
desig = repmat(desig,1,d);
nd = size(desig,1);  %location of design point
n0 = 10;
va = 2; sd = sqrt(va);
span1 = 3; span2 = 1;    % test interval [1,10]
scen = 2.5; sva = 1; ssd = sqrt(sva); ssd2 = ssd;    % x^* distribution
Lm = 100;    % number of macro-replications for random covariates point

m = [5,11,23,49,103,220,470,1000]; nm = length(m);    % number of randomly selected covariate
%%% PFS
testa = 100000;    % number of test covariate point for PFS
Del = 0.2;  %0.2
%%%%%%%
options = optimoptions('fmincon','OptimalityTolerance',0.01,...
    'StepTolerance',0.01,'Display','off');
optseed = 10;  
%%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
[ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0 );

%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Tmu0=2.5; Trusd0 = 3;
Truva = 0.8;
[ TTPmse1,TTPpfs1,TTPmse2,TTPpfs2,TTPmse3,TTPpfs3,TTPmse4,TTPpfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu0,Trusd0 );

%%%%%%%%%%%%%%%%%%%% Sample distribution: Norm %%%%%%%%%%%%%%%%%%%%%%
[ NPmse1,NPpfs1,NPmse2,NPpfs2,NPmse3,NPpfs3,NPmse4,NPpfs4 ] ...
    = NorS( fun,m,desig,scen,ssd,Lm,testa,Del,sd,n0,ssd2 );

%%%%%%%%%%%%%%%%%%%% Adaptive %%%%%%%%%%%%%%%%%%%%%%
m1 = [5,11,23,49,103,220];
Tmu = 2.5; Trusd = 0.75;   
pd = makedist('Normal','mu',Tmu,'sigma',Trusd); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);
Etestx = random(Tpd,testa,d); 
Etesty = fun(desig,Etestx);    % exact testfun value
Eargm = min(Etesty,[],2);    % true optimal design
[ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnS1( fun,m1,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu,Trusd,Etestx,Etesty );

kernind = 1;
[ EPmse1,EPpfs1,Econshortrec1,TT1 ] ...
    = AdaS10D( fun,m1,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 2;
[ EPmse2,EPpfs2,Econshortrec2,TT2 ] ...
    = AdaS10D( fun,m1,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 3;
[ EPmse3,EPpfs3,Econshortrec3,TT3 ] ...
    = AdaS10D( fun,m1,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 4;
[ EPmse4,EPpfs4,Econshortrec4,TT4 ] ...
    = AdaS10D( fun,m1,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);

save test10D;

figure
subplot(1,2,1)
plot(m1,TPmse1,'*-')
hold on;
plot(m1,TPmse2,'o-')
plot(m1,TPmse3,'+-')
plot(m1,TPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Truncated Normal','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m1, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;



subplot(1,2,2)
plot(m1,EPmse1,'*-')
hold on;
plot(m1,EPmse2,'o-')
plot(m1,EPmse3,'+-')
plot(m1,EPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m1, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')





figure
subplot(1,2,1)
plot(m1,TPpfs1,'*-')
hold on;
plot(m1,TPpfs2,'o-')
plot(m1,TPpfs3,'+-')
plot(m1,TPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Truncated Normal','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m1, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;


subplot(1,2,2)
plot(m1,EPpfs1,'*-')
hold on;
plot(m1,EPpfs2,'o-')
plot(m1,EPpfs3,'+-')
plot(m1,EPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m1, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')



figure
subplot(1,3,1)
plot(m,UPmse1,'*-')
hold on;
plot(m,UPmse2,'o-')
plot(m,UPmse3,'+-')
plot(m,UPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Uniform Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;

subplot(1,3,2)
plot(m,TTPmse1,'*-')
hold on;
plot(m,TTPmse2,'o-')
plot(m,TTPmse3,'+-')
plot(m,TTPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Truncated Normal','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,3,3)
plot(m,NPmse1,'*-')
hold on;
plot(m,NPmse2,'o-')
plot(m,NPmse3,'+-')
plot(m,NPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Normal Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')



figure
subplot(1,3,1)
plot(m,UPpfs1,'*-')
hold on;
plot(m,UPpfs2,'o-')
plot(m,UPpfs3,'+-')
plot(m,UPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Uniform Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;



subplot(1,3,2)
plot(m,TTPpfs1,'*-')
hold on;
plot(m,TTPpfs2,'o-')
plot(m,TTPpfs3,'+-')
plot(m,TTPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Truncated Normal','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,3,3)
plot(m,NPpfs1,'*-')
hold on;
plot(m,NPpfs2,'o-')
plot(m,NPpfs3,'+-')
plot(m,NPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (10D): Normal Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')












