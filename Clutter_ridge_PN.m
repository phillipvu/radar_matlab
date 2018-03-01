% clear all
%close all

rng('shuffle')

%% Initialization

N = 10;
M = 1;
P = 16;
Z = 16;
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

A_ula = generate_dictionary(1.5*(0:N-1)',zt,g,P);
C_ula = generate_clutter(2.5*(0:N-1)',zt,1001,P);


%% Generate covariance matrix

R_ra = var_c*(C_ra*C_ra') + var_n*eye(N*P);
R_ula = var_c*(C_ula*C_ula') + var_n*eye(N*P);

P_ra = [];
P_ula = [];

for i = 1:length(g)^2

P_ra(i) = abs(A_ra(:,i)'*R_ra*A_ra(:,i));
P_ula(i) = abs(A_ula(:,i)'*R_ula*A_ula(:,i));

end

P_ra = reshape(P_ra,length(g),length(g));
P_ula = reshape(P_ula,length(g),length(g));

figure, imagesc(g,g,P_ra)
axis([-1 1 -1 1])
