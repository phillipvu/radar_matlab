% clear all
%close all

rng('shuffle')

%% Initialization

N = 8;
M = 1;
P = 16;
Z = 8;
K = 1;

var_n = 10^-(2);
var_c = var_n*1e3;

L = 0.5*(0:P-1)';
g = -1:1/Z/8:1;

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
C_ula1 = generate_clutter(0.5*(0:N-1)',zt,1001,P);

A_ula2 = generate_dictionary(0.5*(0:2*Z-1)',zt,g,P);
C_ula2 = generate_clutter(0.5*(0:2*Z-1)',zt,1001,P);

%% Generate covariance matrix

R_ra = inv(var_c*(C_ra*C_ra') + var_n*eye(N*P));
R_ula1 = inv(var_c*(C_ula1*C_ula1') + var_n*eye(N*P));
R_ula2 = inv(var_c*(C_ula2*C_ula2') + var_n*eye(P^2));

P_ra = [];
P_ula1 = [];
P_ula2 = [];

for i = 1:length(g)^2

P_ra(i) = 10*log10(abs(A_ra(:,i)'*R_ra*A_ra(:,i)));
P_ula1(i) = 10*log10(abs(A_ula1(:,i)'*R_ula1*A_ula1(:,i)));
P_ula2(i) = 10*log10(abs(A_ula2(:,i)'*R_ula2*A_ula2(:,i)));

end

P_ra = reshape(P_ra,length(g),length(g));
P_ula1 = reshape(P_ula1,length(g),length(g));
P_ula2 = reshape(P_ula2,length(g),length(g));

figure, plot(g,P_ra(ceil(length(g)/2),:),g,P_ula1(ceil(length(g)/2),:),g,P_ula2(ceil(length(g)/2),:))
axis([-0.2 0.2 10 20])


