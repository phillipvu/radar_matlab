% clear all
close all

%% Intialize

N_1 = 16;
N_2 = 8;
M = 1;
P = 16;
Z = 8;

var_n = 10^(-2);
CNR = var_n*1e3;

L = 0.5*(0:P-1)';
z_u = 0.5*(0:N_1-1)';
z_r = Z*sort(rand(N_1,1)); z_r(1) = 0; z_r(end) = Z;

Zt = Z/2; Zr = Z/2;

zt = Zt*sort(rand(M,1)); zt(1) = 0; zt(end) = Zt;

zr = Zr*sort(rand(N_2,1)); zr(1) = 0; zr(end) = Zr;
zr = zr + rand;
g = -1:1/Z/16:1;

%% Initialize steering vectors and Clutter covariance matrix

A_u = [];
A_r = [];
A_m = [];

C_u = [];
C_r = [];
C_m = []; 
b_ = [];
for i = 1:length(g)
    b = exp(1i*2*pi*L*g(i));
    
    for j = 1:length(g)
        a_u = exp(1i*2*pi*z_u*g(j));
        a_r = exp(1i*2*pi*z_r*g(j));
        
        t_m = exp(1i*2*pi*zt*g(j));
        r_m = exp(1i*2*pi*zr*g(j));    
        
        a_m = kron(t_m,r_m);
        
        
        A_u(:,end+1) = kron(b,a_u)/sqrt(N_1*P);
        A_r(:,end+1) = kron(b,a_r)/sqrt(N_1*P);
        A_m(:,end+1) = kron(b,a_m)/sqrt(N_2*M*P);
        
        
        if j == i
            C_u(:,end+1) = A_u(:,end);
            C_r(:,end+1) = A_r(:,end);
            C_m(:,end+1) = A_m(:,end);
        end
    end
end

R_u = CNR*C_u*C_u' + var_n*eye(N_1*P);
R_r = CNR*C_r*C_r' + var_n*eye(N_1*P);
R_m = CNR*C_m*C_m' + var_n*eye(N_2*M*P);

% s_u = 10*log10(svd(R_u));
% s_r = 10*log10(svd(R_r));
% s_m = 10*log10(svd(R_m));
% 
% figure, plot(1:2*Z*P,s_u,1:N_1*P,s_r,1:N_2*M*P,s_m)
% axis([0 70 -15 25])

P_u = [];
P_r = [];
P_m = [];

for i = 1:length(g)^2
    P_u(end+1) = real(A_u(:,i)'*R_u*A_u(:,i));
    P_r(end+1) = real(A_r(:,i)'*R_r*A_r(:,i));
    P_m(end+1) = real(A_m(:,i)'*R_m*A_m(:,i));
end


P_u = reshape(P_u,length(g),length(g));
P_r = reshape(P_r,length(g),length(g));
P_m = reshape(P_m,length(g),length(g));

figure, imagesc(g,g,P_u)
xlabel('Spatial frequency')
ylabel('Normalized Doppler')
axis([-0.5 0.5 -0.5 0.5])

figure, imagesc(g,g,P_r)
xlabel('Spatial frequency')
ylabel('Normalized Doppler')
axis([-0.5 0.5 -0.5 0.5])

figure, imagesc(g,g,P_m)
xlabel('Spatial frequency')
ylabel('Normalized Doppler')
axis([-0.5 0.5 -0.5 0.5])

% P_u = [];
% P_r = [];
% P_m = [];
% 
% for i = 1:length(g)^2
%     P_u(end+1) = abs(A_u(:,i)'*A_u(:,600));
%     P_r(end+1) = abs(A_r(:,i)'*A_r(:,600));
%     P_m(end+1) = abs(A_m(:,i)'*A_m(:,600));
% end
% 
% P_u = reshape(P_u,length(g),length(g));
% P_r = reshape(P_r,length(g),length(g));
% P_m = reshape(P_m,length(g),length(g));
% 
% figure, imagesc(g,g,P_u)
% xlabel('Spatial frequency')
% ylabel('Normalized Doppler')
% % axis([-0.5 0.5 -0.5 0.5])
% 
% figure, imagesc(g,g,P_r)
% xlabel('Spatial frequency')
% ylabel('Normalized Doppler')
% % axis([-0.5 0.5 -0.5 0.5])
% 
% figure, imagesc(g,g,P_m)
% xlabel('Spatial frequency')
% ylabel('Normalized Doppler')
% % axis([-0.5 0.5 -0.5 0.5])

