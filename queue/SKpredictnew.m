function [f, MSE] = SKpredictnew(model,Xpred,Bpred)
% make predictions at prediction points using a stochastic kriging model  
% model = output of SKfit
% Xpred = (K x d) matrix of prediction points
% f = (K x 1) predictions at predictions points
% 
% Exmaples
%      SK_gau  = SKpredict(skriging_model,XK,ones(K,1));
% Based on parameter estimates skriging_model obtained from SKfit.m,
% use SK model to predict the values at prediction points XK with constant
% prediction trend, Bpred = ones(K,1)

% retrieve model parameters from model structure obtained from SKfit
X = model.X;
% minX = model.minX;
% maxX = model.maxX;
[k d] = size(X);
theta = model.theta;
kernind = model.kernind;
Znew = model.Znew;
Sigmainv = model.Sigmainv;
beta = model.beta;
tau = sqrt(model.tausquared);

% simple check for dimensions of Xpred and X
K = size(Xpred,1);     % number of prediction points

% % calculate distance matrix for prediction points
% Xpred = (Xpred - repmat(minX,K,1)) ./ repmat(maxX-minX,K,1);

% calculate correlations between prediction points and design points
Rpred = Kernval(Xpred, X, [theta,tau], kernind);

% calculate responses at prediction points 
f = Bpred*beta + Rpred'*Znew;

if nargout > 1        
    MSE = zeros(K,1);
    B = model.B;
    for i = 1 : K
        r = Rpred(:,i);        
%         MSE(i,1) = tau2 - r'*Sigmainv*r;
        delta = Bpred(i,:)' - B'*Sigmainv*r;    
        MSE(i,1) = tau2 - r'*Sigmainv*r + delta'*(B'*Sigmainv*B)^(-1)*delta;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%