% clear all
close all

rng('shuffle')

%% Initialization

N = 4;
M = 2;
P = 16;
Z = 8;
K = 2;
fau = logspace(-10,0,10);
fa = logspace(-10,-2,10);

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z:1;

S =  [9 -25] + length(g)*9;
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

P_map = abs(A_ra'*R_c*A_ra(:,S(1))).^2./A_';
P_map2 = sort(P_map,'descend');

mu = 10*log10(P_map2(2)/P_map2(1))
P_map = reshape(P_map,17,17);

if mu > PSL && mu < PSL+1
    break
end

end

%% Monte Carlo Simulations

ex = 5e4;

Pd = zeros(length(fa),4);
Pf = Pd; 

A_raS = A_ra(:,S);
    

del = 1/2/Z*rand;
u =  [del del];

f_ = 1/Z + del;
f_ = [-f_ f_];

A_ra2 = [];

for s = 1:K
    b = exp(1i*2*pi*L*f_(s));
    t = exp(1i*2*pi*zt*u(s));
    r = exp(1i*2*pi*zr*u(s));
    A_ra2(:,end+1) = kron(b,kron(t,r))/sqrt(N_);
end    

for i = 1:ex

clc
[i/ex mu]
Pd/(i-1)/K
Pf/(i-1)
     
% Generate secondary data

Z_ra = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_) + randn(size(C_ra,2),2*N_)*1i);

% Generate primary data

y_ra = A_raS*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);
y_ra2 = A_ra2*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);


% Perform detection

[S_AMF2,P2] = CFAR_AMF(y_ra,A_ra,Z_ra,fau);
S_MP = CFAR_MHMP4(y_ra,A_ra,[],fa,Z_ra);
S_MP2 = CFAR_MHMP5(y_ra2,zr,zt,P,g,1/8/Z,[],fa,Z_ra);
S_MP3 = CFAR_MHMP5(y_ra2,zr,zt,P,g,1/8/Z,[5 1 1],fa,Z_ra);

% Record results
   

for f = 1:length(fa)
    
    [a,b] = off_grid([f_',u'],S_MP2(f).AD,sqrt(2)/Z); 
    [a2,b2] = off_grid([f_',u'],S_MP3(f).AD,sqrt(2)/Z);
    
    Pd(f,:) = Pd(f,:) + [length(intersect(S_AMF2(f).S,S)) length(intersect(S_MP(f).S,S)) a a2];
    Pf(f,:) = Pf(f,:) + 1 - [isempty(setdiff(S_AMF2(f).S,S)) isempty(setdiff(S_MP(f).S,S)) (b==0) (b2==0)];  
    
end

% if b>0
%     clc
% end


end

Pd = Pd/ex/K;
Pf = Pf/ex;



% figure, semilogx(Pf(:,1),Pd(:,1),Pf(:,2),Pd(:,2),'-rx',Pf(:,3),Pd(:,3),'-ks',Pf(:,4),Pd(:,4),'-go','LineWidth',2,'MarkerSize',10)
% xlabel('P_F')
% ylabel('P_D')
% legend('ABF - (4\lambda ULA)','ABF - (8\lambda ULA)','ABF - (8\lambda RA','MP-CFAR -  (8\lambda RA)')


% figure, semilogx(Pf(:,1),Pd(:,1),Pf(:,2),Pd(:,2),'-rx',Pf(:,3),Pd(:,3),'-ks',Pf(:,4),Pd(:,4),'-go','LineWidth',2,'MarkerSize',10)
% xlabel('P_F')
% ylabel('P_D')
% legend('ABF - (4\lambda ULA)','ABF - (8\lambda ULA)','ABF - (8\lambda RA','MP-CFAR - (8\lambda RA)', 'MBMP-CFAR - (8\lambda RA)')


