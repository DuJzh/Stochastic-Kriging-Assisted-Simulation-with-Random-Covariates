
clearvars,clc;

% Griewangk's function: testfun2
fun = @testfun2;
desig = 1:10; desig = desig'; d = 10;  % context dim.  desig=1:4;
desig = repmat(desig,1,d);
nd = size(desig,1);  %location of design point
n0 = 10;
va = 2; sd = sqrt(va);
span1 = 3; span2 = 1;    % test interval [1,10]
Lm = 100;    % number of macro-replications for random covariates point

%%% PFS
testa = 100000;    % number of test covariate point for PFS
Del = 0.2;  %0.2
%%%%%%%
options = optimoptions('fmincon','OptimalityTolerance',0.01,...
    'StepTolerance',0.01,'Display','off');
optseed = 10;  
%%%%%%%%%%%%%%%%%%%% Adaptive %%%%%%%%%%%%%%%%%%%%%%
m = [5,11,23,49,103,220];
Etestx = rand(testa,d)*span1+span2;  
Etesty = fun(desig,Etestx);    % exact testfun value
Eargm = min(Etesty,[],2);    % true optimal design
[ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniS1( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Etestx,Etesty );

kernind = 1;
[ EPmse1,EPpfs1,Econshortrec1,TT1 ] ...
    = AdaS10D( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 2;
[ EPmse2,EPpfs2,Econshortrec2,TT2 ] ...
    = AdaS10D( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 3;
[ EPmse3,EPpfs3,Econshortrec3,TT3 ] ...
    = AdaS10D( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 4;
[ EPmse4,EPpfs4,Econshortrec4,TT4 ] ...
    = AdaS10D( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);

save test10Dadap_unif;

figure
subplot(1,2,1)
plot(m,UPmse1,'*-')
hold on;
plot(m,UPmse2,'o-')
plot(m,UPmse3,'+-')
plot(m,UPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;



subplot(1,2,2)
plot(m,EPmse1,'*-')
hold on;
plot(m,EPmse2,'o-')
plot(m,EPmse3,'+-')
plot(m,EPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')





figure
subplot(1,2,1)
plot(m,UPpfs1,'*-')
hold on;
plot(m,UPpfs2,'o-')
plot(m,UPpfs3,'+-')
plot(m,UPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
% set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;


subplot(1,2,2)
plot(m,EPpfs1,'*-')
hold on;
plot(m,EPpfs2,'o-')
plot(m,EPpfs3,'+-')
plot(m,EPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewangk (10D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')








