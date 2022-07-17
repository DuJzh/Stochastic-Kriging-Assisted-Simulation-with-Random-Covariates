% Purpose: Use stochastic kriging to fit the response surface of expected 
%          waiting time in M/M/1 queue 
% Variable Definition:
%       X - Design points
%       Y - Simulation outputs at design points
%       Vhat - Intrinsic variance at design points
%       k - number of design points
%       K - number of prediction points

clc; clear all; close all;

% === >>> Generate evenly distributed design and prediction points:
maxx = 2; minx = 1.1;   % range of utilization
arrival_rate = 1;       % fixed arrival rate
K = 1000;               % number of prediction points 
k = 10;                 % number of design points
runlength = 3000;       % runlength at each design point
C = 960;                % total computation budget
X = (minx:((maxx-minx)/(k-1)):maxx)';               % design points
truek = arrival_rate./(X.*(X-arrival_rate));        % analytic values at design points
XK = (minx:((maxx-minx)/(K-1)):maxx)';              % prediction points
true = arrival_rate./(XK .* (XK - arrival_rate));   % analytic values at prediction points
 
% === >>> Obtain outputs at design points:
% Effort allocation is proportional to standard deviation
rho = 1./X;
ratio = sqrt(4*rho./(1-rho).^4);
n = ceil(C*ratio/sum(ratio));   % replications at each design point
q = 0;                          % degree of polynomial to fit in regression part(default)
B = repmat(X,[1 q+1]).^repmat(0:q,[k 1]);       % basis function matrix at design points
BK = repmat(XK,[1 q+1]).^repmat(0:q,[K 1]);     % basis function matrix at prediction points
[Y Vhat] = MM1sim(X,arrival_rate,n,runlength,'stationary');  % simulate M/M/1 queue

% === >>> SK fit and predict:
% stochastic kriging with constant trend (q=0), gammaP=2(use gauss correlation function)
skriging_model_2 = SKfit(X,Y,B,Vhat,2);
SK_gau  = SKpredict(skriging_model_2,XK,BK);
    
% === >>> plot SK fitted surface:
titlefontsize = 14;
fontsize = 12;
linewidth = 2;

figure;
plot(XK,SK_gau,'r-','LineWidth',linewidth);
hold on;
plot(XK,true,'k','LineWidth',linewidth);
myleg = legend('stochastic kriging (gauss)',...
        'true values', ...
        'Location','SouthEast');
plot(XK,zeros(K,1),'k-');
ylabel('average expected waiting time');
xlabel('service rate');
hold off;

