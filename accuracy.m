% clear all
% close all

rng('shuffle')

%% Initialization

N = 4;
M = 2;
P = 16;
Z = 8;
K = 1;

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z:1;

mu = inf;

PSL = 10*log10(log(length(g))/N);

N_ = N*(P)*M;
n = 2.^(0:1:4);

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

A_ula1 = generate_dictionary(0.5*(0:N-1)',zt,g,P+1);
C_ula1 = generate_clutter(0.5*(0:N-1)',zt,201,P+1);

A_ula2 = generate_dictionary(0.5*(0:2*Z)',zt,g,P+1);
C_ula2 = generate_clutter(0.5*(0:2*Z)',zt,201,P+1);

R_c = inv(var_c(end)*C_ra*C_ra' + var_n(end)*eye(N_));
A_ = sum(abs(R_c*A_ra).^2,1);

[P_,I_] = sort(A_);
I_ = I_(1:19);

P_map = abs(A_ra'*R_c*A_ra(:,S(1))).^2./A_';
P_map2 = sort(P_map,'descend');

mu = 10*log10(P_map2(2)/P_map2(1))
P_map = reshape(P_map,17,17);

if mu > PSL-1 && mu < PSL
    break
end

end



%% Monte Carlo Simulations

ex = 1e4;
acc = zeros(length(n),ex);

u_ = [];
v_ = [];
    
temp = zeros(length(n),1);    
for i = 1:ex

clc 
i/ex
sum(acc,2)/(i-1)


u =  0 + 1/Z*rand;
f_ = 1/Z + 1/Z*rand;
A_ra2 = [];

% u_(end+1) = u;
% v_(end+1) = f_;

for s = 1:length(S)
    b = exp(1i*2*pi*L*f_(1));
    t = exp(1i*2*pi*zt*u(s));
    r = exp(1i*2*pi*zr*u(s));
    A_ra2(:,end+1) = kron(b,kron(t,r))/sqrt(N_);
end    
    
    
% Generate secondary data

Z_ra = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_) + randn(size(C_ra,2),2*N_)*1i);
 
% Generate primary data

y_ra = A_ra2*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);

% Perform detection
[f_, u]
S_MP = CFAR_MHMP7(y_ra,zr,zt,P,g,1./n/Z,[],fa,Z_ra);

% Record results

for j = 1:length(n)
    acc(j,i) = acc_err([f_ u],S_MP(j).AD); 
end
      

end

