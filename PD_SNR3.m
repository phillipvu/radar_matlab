% clear all
close all

rng('shuffle')
%% Initialization

N = 8;
M = 1;
P = 16;
Z = 8;
K = 1;
fa = 1e-3;

var_n = 10.^(-1:-0.1:-2);
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

S_off = [0 -1 -2];
Sx = setdiff(-9:6,S_off);
I = length(g)*8 + 10;

ex = 2e4;

g = (fa^(-1/(2*N_-N_+2))-1)*(2*N_+1)/(2*N_-N_+1);
Pd = zeros(length(var_n),1);

for v = 1:length(var_n)

a = zeros(1,ex);
b = a;
    
parfor i = 1:ex

y_ = sqrt(var_n(v)/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c(v)/2)*C*(randn(size(C,2),2*N_) + randn(size(C,2),2*N_)*1i);

if K > 0
    x = randperm(length(Sx),1);
    S = I+Sx(x);
    
    y = A(:,[I S])*exp(1i*2*pi*rand(K+1,1)) + sqrt(var_n(v)/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c(v)/2)*C*(randn(size(C,2),1) + randn(size(C,2),1)*1i);
    a(i) = AMF_CFAR(y,A(:,S),A(:,I),y_);
else
    y = A(:,I)*exp(1i*2*pi*rand(1,1)) + sqrt(var_n(v)/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c(v)/2)*C*(randn(size(C,2),1) + randn(size(C,2),1)*1i);
    a(i) = AMF_CFAR(y,A(:,I),[],y_);
end

end

Pd(v) = sum(a>g)/ex;

end

figure, plot(-10*log10(var_n),Pd,'LineWidth',2)
xlabel('SNR(dB)')
ylabel('P_D')




