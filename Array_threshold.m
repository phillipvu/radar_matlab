% clear all
close all

rng('shuffle')

%% Initialization

N = 4;
M = 2;
P = 16;
Z = 8;
K = 1;
fa = 1e-3;
% fau = logspace(-8,-1,8);

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/16/Z:1;
g_ = length(g);


mu = inf;

PSL = 10*log10(log(length(g))/N);

N_ = N*P*M;
N_2 = P*M*length(0.5*(0:2*Z));

%% Generate array

    
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

A_ula2 = generate_dictionary(0.5*(0:2*Z)',zt,g,P);
C_ula2 = generate_clutter(0.5*(0:2*Z)',zt,201,P);


    
R_cu = inv(var_c(end)*C_ula2*C_ula2' + var_n(end)*eye(N_2));
L_ = 2*N_;
L_2 = 2*N_2;
fa_ = 1-(1-fa).^(1/size(A_ra,2));

xt = (L_+1)/(L_-N_+1)*(fa_.^(-1/(L_-N_+2))-1);
xt2 = (L_2+1)/(L_2-N_2+1)*(fa_.^(-1/(L_2-N_+2))-1);



S =  randi(g_^2,1);
A_raS = A_ra(:,S);
A_ula1S = A_ula1(:,S);
A_ula2S = A_ula2(:,S);   

Z_ula2 = sqrt(var_n/2)*(randn(N_2,2*N_2) + randn(N_2,2*N_2)*1i) + sqrt(var_c/2)*C_ula2*(randn(size(C_ula2,2),2*N_2) + randn(size(C_ula2,2),2*N_2)*1i);
y_ula2 = A_ula2S*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_2,1) + randn(N_2,1)*1i) + sqrt(var_c/2)*C_ula2*(randn(size(C_ula2,2),1) + randn(size(C_ula2,2),1)*1i);
[S_AMF2,P2] = CFAR_AMF(y_ula2,A_ula2,Z_ula2,fa);
P2 = reshape(P2,g_,g_);

figure, plot(g,P2(:,floor(S/g_)),g,xt2*ones(length(g),1),'--r')
xlabel('Spatial frequency')
ylabel('ABF output') 

figure, imagesc(g,g,P2)
xlabel('Normalized Doppler')
ylabel('Spatial frequency')




Z_ra = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),2*N_) + randn(size(C_ra,2),2*N_)*1i);
y_ra = A_raS*exp(1i*2*pi*rand(K,1)) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C_ra*(randn(size(C_ra,2),1) + randn(size(C_ra,2),1)*1i);
[S_AMF3,P3] = CFAR_AMF(y_ra,A_ra,Z_ra,fa*10);

P3 = reshape(P3,g_,g_);

figure, plot(g,P3(:,floor(S/g_)),g,xt*ones(length(g),1),'--r')
xlabel('Spatial frequency')
ylabel('Output of ABF') 

figure, imagesc(g,g,P3)
xlabel('Normalized Doppler')
ylabel('Spatial frequency')



