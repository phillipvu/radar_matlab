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

while 1
zr = sort(Z*rand(N,1)); 
zr = zr-min(zr);
if min(diff(zr)) > 0.5
    break
end
end

g = -1:1/Z:1;
N_ = N*P;

%% Initialize steering vectors and Clutter covariance matrix

A = generate_dictionary(zr,0,g,P);
C = generate_clutter(zr,0,401,P);

R_c = inv(var_c(end)*C*C' + var_n(end)*eye(N_));

A_ = sum(abs(R_c*A).^2,1);

[P_,I_] = sort(A_);
I_ = I_(1:19);

ex = 1e4;

PF_1 = zeros(length(var_n),1);
PF_2 = PF_1;
% I = 0;
T_off = zeros(length(var_n),ex);

while 1
   I = randperm(length(g)^2,1);

   if isempty(intersect(I,I_)) == 1
       break
   end
end


P_map = abs(A'*R_c*A(:,I)).^2./A_';
P_map2 = sort(P_map,'descend');


10*log10(P_map2(1)/P_map2(2))
P_map = reshape(P_map,17,17);


for v = 7:7
   a = zeros(1,ex);
   b = a;
    
   for i = 1:ex 
    clc
    v
    i/ex
    [sum(a) sum(b)]/ex
    [PF_1 PF_2]



    y_ = sqrt(var_n(v)/2)*(randn(N_,2*N_) + randn(N_,2*N_)*1i) + sqrt(var_c(v)/2)*C*(randn(size(C,2),2*N_) + randn(size(C,2),2*N_)*1i);
    y = A(:,I)*exp(1i*2*pi*rand(K,1)) + sqrt(var_n(v)/2)*(randn(N_,1) + randn(N_,1)*1i) + sqrt(var_c(v)/2)*C*(randn(size(C,2),1) + randn(size(C,2),1)*1i);

    [S_AMF,P] = CFAR_AMF(y,A,y_,fa);
    a(i) = 1 - isempty(setdiff(S_AMF.S,I));
    
    S_MP = CFAR_MHMP4(y,A,[],fa,y_);
    b(i) = 1 - isempty(setdiff(S_MP.S,I));
      
    if b(i) > 0 && S_MP.S(end) ~= I
%         clc
        T_off(v,i) = 1;
    end
    
   end
   
   PF_1(v) = sum(a)/ex;
   PF_2(v) = sum(b)/ex;
   
end

figure, semilogy(-10*log10(var_n),PF_1,-10*log10(var_n),PF_2,'-kx',-10*log10(var_n),ones(length(var_n),1)*1e-3,'--r','LineWidth',2,'MarkerSize',10)
xlabel('SNR (dB)')
ylabel('P_{FA}')

