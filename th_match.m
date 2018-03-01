% clear all
close all

rng('shuffle')

%% Initialization

N = 4;
M = 2;
P = 16;
Z = 8;
K = 1;
fa = logspace(-1,-3,10);

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z:1;

S =  [9] + length(g)*9;
mu = inf;

N_ = N*P*M;
PSL = 10*log10(log(length(g))/N/M);

%% Generate array

while 1
    
if M > 1
    Zt = Z/2; Zr = Z/2;
    zt = Zt*sort(rand(M,1)); zt(1) = 0; 
    zr = Zr*sort(rand(N,1)); 
else    
    Zr = Z;
    zt = 0;
    zr = Zr*sort(rand(N,1)); zr(1) = 0;    
end

A_ra = generate_dictionary(zr,zt,g,P);
C_ra = generate_clutter(zr,zt,201,P);

R_c = inv(var_c(end)*C_ra*C_ra' + var_n(end)*eye(N_));
A_ = sum(abs(R_c*A_ra).^2,1);

[P_,I_] = sort(A_);
I_ = I_(1:19);

P_map = abs(A_ra'*R_c*A_ra(:,S(1))).^2;
P_map2 = sort(P_map,'descend');

mu = 10*log10(P_map2(2)/P_map2(1))
P_map = reshape(P_map,17,17);

if mu < PSL-1
    break
end

end

%% Monte Carlo Simulations

ex = 1e4;
x = sqrt(P_map2(1));

Pd = zeros(length(fa),1);
Pf = Pd; 

A_raS = A_ra(:,S);
    

for i = 1:ex

clc
[i/ex mu]
% Pd/(i-1)/K
[Pf/(i-1) fa']
     
% Generate secondary data

Z_ra = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_) + randn(size(C_ra,2),2*N_)*1i);

% Generate primary data

y_ra = A_raS*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);


% Perform detection

% [S_AMF,P2] = STAP_BF(y_ra,A_ra,Z_ra,t);
S_MP2 = MBMP_STAP2(y_ra,A_ra,[],fa,R_c);


% Record results
   

for f = 1:length(fa)
        
    Pd(f,:) = Pd(f,:) + length(intersect(S_MP2(f).S,S));
    Pf(f,:) = Pf(f,:) + 1 - isempty(setdiff(S_MP2(f).S,S));  
    
end

% if b>0
%     clc
% end
    

end

Pd = Pd/ex/K;
Pf = Pf/ex;



