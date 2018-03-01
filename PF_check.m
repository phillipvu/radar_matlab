% clear all
close all

rng('shuffle')
%% Initialization

N = 8;
M = 1;
P = 16;
Z = 8;
K = 1;
fa = logspace(-3,-1,6);

var_n = 1e-2;
var_c = var_n*1e3;

L = 0.5*(0:P-1)';

zt = 0;
zr = Z*rand(N,1); 
zr = zr-min(zr);

g = -1:1/Z:1;
N_ = N*P;

%% Initialize steering vectors and Clutter covariance matrix

A = generate_dictionary(zr,0,g,P);
C = generate_clutter(zr,0,401,P);

I = randperm(length(g)^2,1);
while 1
S = randperm(length(g)^2,1);

if isempty(intersect(S,I)) == 1
    break
end

end

ex = 1e4;
L_ = 2*N_;
g1 = (L_+1)/(L_-N_+1)*(fa.^(-1/(L_-N_+2))-1)';



PF_1 = zeros(length(fa),ex);
PF_2 = PF_1;
PF_3 = PF_1;

% T = zeros(1,ex);
% S_MP = []

for i = 1:ex
    
clc
i/ex
[sum(PF_1,2)/(i-1) sum(PF_2,2)/(i-1) sum(PF_3,2)/(i-1) fa']
    
y_ = sqrt(var_n/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c/2)*C*(randn(size(C,2),2*N_) + randn(size(C,2),2*N_)*1i);

y = sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C*(randn(size(C,2),1) + randn(size(C,2),1)*1i);
y_I = A(:,I)*exp(1i*2*pi*rand) + sqrt(var_n/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c/2)*C*(randn(size(C,2),1) + randn(size(C,2),1)*1i);

Z = inv(y_*y_');
T = abs(A(:,S)'*Z*y)^2/real(A(:,S)'*Z*A(:,S));
PF_1(:,i) = T>=g1;

S_MP = CFAR_MHMP4(y,A,[],fa,y_);
S_MP1 = CFAR_MHMP4(y_I,A,[],fa,y_);

for f = 1:length(fa)
    
    if isempty(S_MP(f).S) == 0
        PF_2(f,i) = 1;
    end    
    
    if isempty(setdiff(S_MP1(f).S,I)) == 0
        PF_3(f,i) = 1;
    end
    
end


end

x = 1-(1-fa).^(1/size(A,2));
g = (x.^(-1./(L_-N_+2))-1).*(L_+1)./(L_-N_+1);

PF1 = sum(PF_1,2)/ex;
PF2 = sum(PF_2,2)/ex;
PF3 = sum(PF_3,2)/ex;

figure, semilogy(g,PF2,'-bx',g,PF3,'-ks',g,fa,'--r',g1,fa,':m','LineWidth',2,'MarkerSize',10)
xlabel('Detection threshold')
ylabel('P_{FA}')
legend('MP CFAR - 0 targets','MP-CFAR - 1 target','Analytical MP-CFAR','Analytical ABF')



