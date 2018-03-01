function [S_,P] = CFAR_AMF(y,A,Z,t)

%% Inputs 

% y - Received measurement vector
% A - Dictionary vector
% Z - Secondary data
% t - Desired false alarm probability

%% Outputs

% S - Resolution cells with detected targets
% P - Beampattern

%% Initialization

R_inv = sqrtm(inv(Z*Z'));
A_ = R_inv*A;
normA = sum(abs(A_).^2,1);
z = R_inv*y;

[N,L] = size(Z);
x = 1-(1-t).^(1/size(A,2));
g = (x.^(-1/(L-N+2))-1)*(L+1)/(L-N+1);

%% Perform detection

P = zeros(size(A,2),1);
for i = 1:size(A,2)
   P(i) = abs(A_(:,i)'*z)^2/normA(i);    
end

S_ = [];

for i = 1:length(t)
S_(end+1).S = find(P >= g(i));
S_(end).t = t(i);
end

% S_(1) = [];

end