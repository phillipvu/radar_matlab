function [S_,P] = STAP_BF(y,A,R_inv,fa)

%% Inputs 

% y - Received measurement vector
% A - Dictionary vector
% R_inv - Inverse of the covariance matrix
% fa - Desired false alarm probability

%% Outputs

% S - Resolution cells with detected targets
% P - Beampattern

%% Initialization

N = size(y,1);
t = chi2inv(1-fa,2);

%% Perform detection

P = zeros(size(A,2),1);

for i = 1:size(A,2)
   P(i) = abs(A(:,i)'*R_inv*y)^2/abs(A(:,i)'*R_inv*A(:,i));    
end

S_ = [];

for i = 1:length(t)
S_(end+1).S = find(P >= t(i));
S_(end).t = t(i);
end

% S_(1) = [];

end