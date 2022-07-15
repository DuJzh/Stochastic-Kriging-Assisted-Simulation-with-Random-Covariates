function [f, beta, A] = neglogPL(parms,k,X,B,V,Y,kernind)

% negative of log profile likelihood  
% This is a function of theta and tau2.
% we have profiled out over beta
% parms = [tausq;theta] 
% k = number of design points
% d = dimension of space
% X = design points
% B = basis functions k x d
% V = intrinsic covariance matrix
% Y = simulation output
    
% sum of extrinsic and intrinsic covariances
Sigma  =  Kernval( X,X,parms,kernind ) + V;
% [U,pd] = chol(Sigma);

%%%%%%%%
A = Sigma;
eigen_A = eig(A);
% f = 0.5*k*log(2*pi) + 0.5*sum(log(eigen_A)) + 0.5*Y'/A*Y;
Ainv = A^(-1);
BAinv = B' * Ainv;
beta = (BAinv*B) \ (BAinv*Y);
f = 0.5*k*log(2*pi) + 0.5*sum(log(eigen_A)) + 0.5*(Y-B*beta)'*Ainv*(Y-B*beta);

% BAinv = B' / A;
% beta = (BAinv*B) \ (BAinv*Y);
% f = 0.5*k*log(2*pi) + 0.5*sum(log(eigen_A)) + 0.5*(Y-B*beta)'/A*(Y-B*beta);
% %%%%%%%%
% 
% % invert it via Cholesky factorization
% L = U';
% Linv = inv(L);
% Sinv = Linv'*Linv;
% 
% % the optimal beta given theta and tau2
% beta = inv(B'*Sinv*B) \ (B'*(Sinv*Y)); 
% Z = L\(Y-B*beta);
% 
% % negative log likelihood
% f = (log(det(L)) + 0.5*Z'*Z + 0.5*k*log(2*pi));

