% clear all
close all

N_1 = 8;
N_2 = 3;
M = 3;
P = 16;
Z = 8;

var_n = 10^(-2);
CNR = var_n*1e3;
    
L = 0.5*(0:P-1)';
z_u = 0.5*(0:N_1-1)';
z_ul = 0.5*(0:2*Z-1)';
z_r = Z*sort(rand(N_1,1)); z_r(1) = 0; z_r(end) = Z;

Zt = 4; Zr = 4;
zt = Zt*sort(rand(M,1)); zt(1) = 0; zt(end) = Zt;
zr = Zr*sort(rand(N_2,1)); zr(1) = 0; zr(end) = Zr;

g = -1:1/Z/16:1;

%% Initialize steering vectors and Clutter covariance matrix

A_u = [];
A_ul = [];
A_r = [];
A_m = [];

C_u = [];
C_ul = [];
C_r = [];
C_m = []; 

for i = 1:length(g)
    b = exp(1i*2*pi*L*g(i));
    for j = 1:length(g)
        a_u = exp(1i*2*pi*z_u*g(j));
        a_ul = exp(1i*2*pi*z_ul*g(j));
        a_r = exp(1i*2*pi*z_r*g(j));
        
        t_m = exp(1i*2*pi*zt*g(j));
        r_m = exp(1i*2*pi*zr*g(j));    
        
        a_m = kron(t_m,r_m);
        
        
        A_u(:,end+1) = kron(b,a_u)/sqrt(N_1*P);
        A_ul(:,end+1) = kron(b,a_ul)/sqrt(2*Z*P);
        A_r(:,end+1) = kron(b,a_r)/sqrt(N_1*P);
        A_m(:,end+1) = kron(b,a_m)/sqrt(N_2*M*P);
        
        
        if j == i
            C_u(:,end+1) = A_u(:,end);
            C_ul(:,end+1) = A_ul(:,end);
            C_r(:,end+1) = A_r(:,end);
            C_m(:,end+1) = A_m(:,end);
        end
    end
end

R_u = CNR*C_u*C_u' + var_n*eye(N_1*P);
R_ul = CNR*C_ul*C_ul' + var_n*eye(2*Z*P);
R_r = CNR*C_r*C_r' + var_n*eye(N_1*P);
R_m = CNR*C_m*C_m' + var_n*eye(N_2*M*P);

s_u = 10*log10(svd(R_u));
s_ul = 10*log10(svd(R_ul));
s_r = 10*log10(svd(R_r));
s_m = 10*log10(svd(R_m));

figure, plot(1:N_1*P,s_u,  1:P*2*Z,s_ul,'-r',  1:5:N_1*P,s_r(1:5:end),'--kx',  1:5:N_2*M*P,s_m(1:5:end),'--ms', 'LineWidth',2,'MarkerSize',10)
axis([0 70 -20 25])
xlabel('Eigenvalue index')
ylabel('Eigenvalue(dB)')
legend('ULA - 4\lambda', 'ULA - 8\lambda', 'RA - 8\lambda', 'MIMO-RA - 8\lambda')

%% SINR

P_u = [];
P_ul = [];
P_r = [];
P_m = [];

a_u = ones(N_1,1);
a_ul = ones(2*Z,1);
a_r = ones(N_1,1);
a_m = ones(N_2*M,1);

R_ui = inv(R_u);
R_uli = inv(R_ul);
R_ri = inv(R_r);
R_mi = inv(R_m);

for i = 1:length(g)
    b = exp(1i*2*pi*L*g(i));
    
    a_u = kron(b,ones(N_1,1))/sqrt(N_1*P);
    a_ul = kron(b,ones(2*Z,1))/sqrt(2*Z*P);
    a_r = kron(b,ones(N_1,1))/sqrt(N_1*P);
    a_m = kron(b,ones(N_2*M,1))/sqrt(N_2*M*P);
    
    P_u(end+1) = real(a_u'*R_ui*a_u);
    P_ul(end+1) = real(a_ul'*R_uli*a_ul);
    P_r(end+1) = real(a_r'*R_ri*a_r);
    P_m(end+1) = real(a_m'*R_mi*a_m);

end

figure, plot(g,10*log10(P_u), g,10*log10(P_ul),'-r', g(1:2:end),10*log10(P_r(1:2:end)),'--kx', g(1:2:end),10*log10(P_m(1:2:end)),'--ms', 'LineWidth',2,'MarkerSize',10)
axis([-0.1 0.1 0 20])
xlabel('Normalized Doppler')
ylabel('SINR(dB)')
legend('ULA - 4\lambda', 'ULA - 8\lambda', 'RA - 8\lambda', 'MIMO-RA - 8\lambda')

figure, plot(g,10*log10(P_u), g,10*log10(P_ul),'-r', g(1:4:end),10*log10(P_r(1:4:end)),'--kx', g(1:4:end),10*log10(P_m(1:4:end)),'--ms', 'LineWidth',2,'MarkerSize',10)
axis([0 1 0 22])
xlabel('Normalized Doppler')
ylabel('SINR(dB)')
legend('ULA - 4\lambda', 'ULA - 8\lambda', 'RA - 8\lambda', 'MIMO-RA - 8\lambda')
