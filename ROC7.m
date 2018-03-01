% clear all
% close all

rng('shuffle')

%% Initialization

N = 8;
M = 1;
P = 16;
Z = 8;
K = 1;
fa = logspace(-12,-10,8);
% fau = logspace(-8,-1,8);

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z:1;

S =  [9] + length(g)*9;
mu = inf;

PSL = 10*log10(log(length(g))/N);

N_ = N*P*M;
N_2 = P*M*length(0.5*(0:2*Z));

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
C_ra = generate_clutter(zr,zt,1001,P);

A_ula1 = generate_dictionary(0.5*(0:N-1)',zt,g,P);
C_ula1 = generate_clutter(0.5*(0:N-1)',zt,201,P);

A_ula2 = generate_dictionary(0.5*(0:2*Z)',zt,g,P);
C_ula2 = generate_clutter(0.5*(0:2*Z)',zt,201,P);

R_c = inv(var_c(end)*C_ra*C_ra' + var_n(end)*eye(N_));
A_ = sum(abs(sqrtm(R_c)*A_ra).^2,1);

[P_,I_] = sort(A_);
I_ = I_(1:19);

P_map = abs(A_ra'*R_c*A_ra(:,S(1))).^2./A_';
P_map2 = sort(P_map,'descend');

mu = 10*log10(P_map2(2)/P_map2(1))
P_map = reshape(P_map,17,17);

% if mu < PSL+1 && mu > PSL-1
    break
% end

end

%% Monte Carlo Simulations

ex = 2e4;

Pd = zeros(length(fa),4);
Pf = Pd; 


A_raS = A_ra(:,S);
A_ula1S = A_ula1(:,S);
A_ula2S = A_ula2(:,S);   
    
R_cu = inv(var_c(end)*C_ula2*C_ula2' + var_n(end)*eye(N_2));
L_ = 2*N_;
L_2 = 2*N_2;
fa_ = 1-(1-fa).^(1/size(A_ra,2));

xt = (L_+1)/(L_-N_+1)*(fa_.^(-1/(L_-N_+2))-1);
xt2 = (L_2+1)/(L_2-N_2+1)*(fa_.^(-1/(L_2-N_+2))-1);


Pdu = AMF_PD5(A_ula2(:,S),A_ula2(:,S),L_2,R_cu,xt2);
Pdr = AMF_PD5(A_ra(:,S),A_ra(:,S),L_,R_c,xt);


for i = 1:ex

clc
[i/ex mu]
[Pd(:,[2 4])/(i-1)/K]
[Pf(:,[2 4])/(i-1)]
    
% Generate secondary data

% Z_ra = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_) + randn(size(C_ra,2),2*N_)*1i);
% Z_ula1 = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ula1*(randn(size(C_ula1,2),2*N_) + randn(size(C_ula1,2),2*N_)*1i);
Z_ula2 = sqrt(var_n/2)*(randn(N_2,2*N_2) + randn(N_2,2*N_2)*1i) + sqrt(var_c/2)*C_ula2*(randn(size(C_ula2,2),2*N_2) + randn(size(C_ula2,2),2*N_2)*1i);

% Generate primary data

% y_ra = A_raS*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);
% y_ula1 = A_ula1S*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ula1*(randn(size(C_ula1,2),1) + randn(size(C_ula1,2),1)*1i);
y_ula2 = sqrt(2)*A_ula2S*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_2,1) + randn(N_2,1)*1i) + sqrt(var_c/2)*C_ula2*(randn(size(C_ula2,2),1) + randn(size(C_ula2,2),1)*1i);

% Perform detection

% [S_AMF1,P1] = CFAR_AMF(y_ula1,A_ula1,Z_ula1,fa*10);
[S_AMF2,P2] = CFAR_AMF(y_ula2,A_ula2,Z_ula2,fa);
% S_MP = CFAR_MHMP4(y_ula2,A_ula2,[],fa,Z_ula2);
% [S_AMF3,P3] = CFAR_AMF(y_ra,A_ra,Z_ra,fa*10);
% S_MP = CFAR_MHMP4(y_ra,A_ra,[],fa,Z_ra);
% S_MBMP = CFAR_MHMP4(y_ra,A_ra,[3 3 1],fa,Z_ra);

% Record results

for f = 1:length(fa)
    
    Pd(f,:) = Pd(f,:) + [length(intersect(S_AMF1(f).S,S)) length(intersect(S_AMF2(f).S,S)) length(intersect(S_AMF3(f).S,S)) length(intersect(S_MP(f).S,S))];
    Pf(f,:) = Pf(f,:) + 1 - [isempty(setdiff(S_AMF1(f).S,S)) isempty(setdiff(S_AMF2(f).S,S)) isempty(setdiff(S_AMF3(f).S,S)) isempty(setdiff(S_MP(f).S,S))]; 
        
end

% if isempty(setdiff(S_AMF2(f).S,S)) == 0 
%     clc
% end


end

Pd = Pd/ex/K;
Pf = Pf/ex;



% figure, semilogx(Pf(:,1),Pd(:,1),Pf(:,2),Pd(:,2),'-rx',Pf(:,3),Pd(:,3),'-ks',Pf(:,4),Pd(:,4),'-go','LineWidth',2,'MarkerSize',10)
% xlabel('P_F')
% ylabel('P_D')
% legend('ABF - (4\lambda ULA)','ABF - (8\lambda ULA)','ABF - (8\lambda RA','MP-CFAR -  (8\lambda RA)')

    
figure, semilogx(Pf(:,1),Pd(:,1),Pf(:,3),Pd(:,3),'-ks',Pf(:,2),Pd(:,2),'-rx',Pf(:,4),Pd(:,4),'-go','LineWidth',2,'MarkerSize',10)
xlabel('P_F')
ylabel('P_D')
legend('ABF - (4\lambda ULA)','ABF - (8\lambda RA)','ABF - (8\lambda ULA)','MP-CFAR', 'MBMP-CFAR')


