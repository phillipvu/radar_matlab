% clear all
close all

rng('shuffle')

%% Initialization

N = 3:7;
M = 2;
P = 16;
Z = 8;
K = 2;
fa = 1e-5;

var_n = 10^-(1.8);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z:1;

S =  [-25 9] + length(g)*9;
mu = inf;

N_ = N*P*M;



%% Monte Carlo Simulations

ex = 5e4;

Pd = zeros(length(N),2);
Pf = Pd; 

A_raS = A_ra(:,S);

Zt = Z/2; Zr = Z/2;    
for n = length(N):-1:1

for i = 1:ex

clc
[n/length(N) i/ex]
Pd/(i-1)/K
    
zt = Zt*sort(rand(M,1)); zt(1) = 0; 
zr = Zr*sort(rand(N(n),1)); 

A_ra = generate_dictionary(zr,zt,g,P);
C_ra = generate_clutter(zr,zt,401,P);


% Generate secondary data

Z_ra = sqrt(var_n/2)*(randn(N_(n),2*N_(n)) + randn(N_(n),2*N_(n))*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_(n)) + randn(size(C_ra,2),2*N_(n))*1i);

% Generate primary data

y_ra = A_raS*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_(n),1) + randn(N_(n),1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);


% Perform detection

S_MP = MBMP_STAP(y_ra,A_ra,ones(10,1),fa,Z_ra);
S_MP2 = MBMP_STAP(y_ra,A_ra,[5 ones(1,9)],fa,Z_ra);

% Record results
           
Pd(n,:) = Pd(n,:) + [length(intersect(S_MP.S,S)) length(intersect(S_MP2.S,S))];
    
end

end

Pd = Pd/ex/K;



