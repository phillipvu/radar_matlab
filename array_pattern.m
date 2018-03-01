% clear all
%close all 

%% Initialization

N = 9;
Z = 9;
g = linspace(-1,1,1003);

zr = Z*rand(N,1); zr = zr-min(zr);
zu1 = 0.5*(0:N-1)';
zu2 = 0.5*(0:2*Z);


A_r = generate_dictionary1(zr,0,g,0);
A_u1 = generate_dictionary1(zu1,0,g,0);
A_u2 = generate_dictionary1(zu2,0,g,0);

s1 = ones(N,1)/sqrt(N);
s2 = ones(2*Z+1,1)/sqrt(2*Z+1);

%% Random array statistics

SL_m = -10*log10(N);
SL_p = 10*log10(log(2*Z)/N);

%% Random array pattern

P_r = 10*log10(abs(A_r'*s1).^2);
% 
% % figure, plot(g,P_r,':r')
% figure, plot(g,P_r,g,ones(length(g),1)*SL_m,'--r',g,ones(length(g),1)*SL_p,'--k','LineWidth',2)
% axis([-1 1 -15 0])
% xlabel('Spatial frequency')
% ylabel('Beampattern (dB)')
% legend('RA beampattern','Mean SL level','Mean PSL level')

%% Uniform linear array beampattern

P_u1 = 10*log10(abs(A_u1'*s1).^2);
P_u2 = 10*log10(abs(A_u2'*s2).^2);

hold on; plot(g,P_u1);
%plot(g,P_u1,g,P_u2,'-r')
 %figure, plot(g,P_u1,g,P_u2,'-r',g,P_r,':k',g,ones(length(g),1)*SL_m,'--r',g,ones(length(g),1)*SL_p,'--g','LineWidth',2)   
 axis([-1 1 -15 0])
 xlabel('Spatial frequency')
 ylabel('Beampattern (dB)')
%legend('ULA - 4\lambda','ULA - 8\lambda','RA - 8\lambda','Average SL level','Average PSL level')
