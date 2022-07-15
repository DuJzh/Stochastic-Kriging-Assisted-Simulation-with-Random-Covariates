function model = SKfitnew(X, Y, B, Vhat, kernind)
% fit a stochastic kriging model to simulation output
% X - design points (k x d) 
% Y - (k x 1) vector of simulation output at each design point
%     The first column must be a column of ones!
% B - (k x b) matrix of basis functions at each design point
% Vhat - (k x 1) vector of simulation output variances
% Types of correlation function used to fit surface:
%    kernind = 1: exponential
%    kernind = 2: squared exponential
%    kernind = 3: matern 3/2
%    kernind = 4: matern 5/2
% 
% Examples
%       skriging_model = SKfit(X,Y,ones(k,1),Vhat,2);
% Use SK model with gauss correlation function to fit data, (X,Y,Vhat)
% X is design points, Y and Vhat are outputs at design points 
% (Y is mean, Vhat is intrinsic variance), non-trend estimate B = ones(k,1)

[k d] = size(X);

% % Normalize data by scaling each dimension from 0 to 1
% minX = min(X);  
% maxX = max(X);
% X = (X - repmat(minX,k,1)) ./ repmat(maxX-minX,k,1);

% diagonal intrinsic variance matrix
V = diag(Vhat);

% initialize parameters theta, tau2 and beta
% inital extrinsic variance = variance of ordinary regression residuals
betahat = (B'*B)\(B'*Y);
tau2_0 = var(Y-B*betahat);

% set initial correlation
ndistX = k*(k-1) / 2;        % max number of non-zero distances
distX = zeros(ndistX, d);    % initialize matrix with distances
temp = 0;
for i = 1 : k-1
    temp = temp(end) + (1 : k-i);
    distX(temp,:) = repmat(X(i,:), k-i, 1) - X(i+1:k,:); 
end
average_distance = mean(sum(abs(distX).^2, 2));
% max_distance = max(sum(abs(distX).^2, 2));
if kernind == 2  % squared exponential
    theta_0 = sqrt(2*average_distance/log(2));
elseif kernind == 1  % exponential
    theta_0 = sqrt(average_distance)/log(2);
elseif kernind == 3  % matern 3/2
    theta_0 = sqrt(3*average_distance)/log(2);
elseif kernind == 4  % matern 5/2
    theta_0 = sqrt(5*average_distance)/log(2);
end

% lower bounds for parameters tau2, theta, and beta included 
% only to satisfy fmincon  
lbtau = sqrt(0.00000001*tau2_0);     % naturally 0 
lbtheta = 0.000001; % naturally 0; increase to avoid numerical trouble
lb = [lbtheta;lbtau];

% maximize profile log-likelihood function (-"logPL")
% subject to lower bounds on tau2 and theta
myopt = optimset('Display','notify','MaxFunEvals',1000000,'MaxIter',500);
parms = fmincon(@(x) neglogPL(x,k,X,B,V,Y,kernind),...
        [theta_0;sqrt(tau2_0)],[],[],[],[],lb,[],[],myopt); 

% record MLEs for tau2 and theta 
tau2hat = (parms(2))^2;
thetahat = parms(1);

% MLE of beta is known in closed form given values for tau2, theta
% and is computed below
% calculate estimates of the correlation and covariance matrices

% issue warnings related to constraints 
warningtolerance = 0.00001;
if min(abs(lbtheta - thetahat)) < warningtolerance
    warning('thetahat was very close to artificial lower bound');
end
if abs(lbtau - parms(2)) < warningtolerance
    warning('tau2hat was very close to artificial lower bound');
end

% output MLEs and other things useful in prediction
model.tausquared =  tau2hat;
model.theta = thetahat;
model.X = X;
% model.minX = minX;
% model.maxX = maxX;
model.kernind = kernind;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output some additional parameters, added by SHEN Haihui May 23 2017 %%%
[negLL, betahat, A] = neglogPL(parms,k,X,B,V,Y,kernind);
model.negLL = negLL;
model.beta = betahat;
model.B = B;
Sigmahat  = A;
Lhat = chol(Sigmahat)';
Lhatinv = inv(Lhat);
Sigmahatinv = Lhatinv'*Lhatinv;
model.Sigma = Sigmahat;
model.Sigmainv = Sigmahatinv;
model.Znew = Sigmahatinv*(Y-B*betahat);

% % AIC, BIC
% AIC = -2*(-model.negLL) + 2*(1+length(betahat)+length(thetahat));
% BIC = -2*(-model.negLL) + (1+length(betahat)+length(thetahat))*log(k);
% statistic = [AIC, BIC];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%