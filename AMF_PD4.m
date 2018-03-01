function Pd = AMF_PD4(x,s,I,L,R_inv,t,c_th)

%% Inputs

% x - Complex amplitude of the target
% s - Steering vector of interest
% I - Interference basis
% L - Number of secondary data used
% R_Inv - Inverse of the theoretical covariance matrix
% t - threshold parameter
% c_th - Mismatch factor

%% Output

% Pd - Probability that the test statistic exceeds the threshold 


%% Initialization

N = length(s);
S = (eye(N) - I*inv(I'*R_inv*I)*I'*R_inv)*s;
M = L-N+1;
x_ = linspace(0,1,10001);


SNR_t = abs(x)^2*2*abs(s'*R_inv*s)^2/abs(s'*R_inv*s);


px = abs(LF_pdf3(L,S,c_th,SNR_t));

Pd = zeros(length(t),1);

for i = 1:length(t)    
    Pd(i) = 1 - sum(ncfcdf(M*t(i)*x_,2,2*M,x_.*SNR_t*c_th).*px)/length(x_);        
end


end

function p = LF_pdf3(L,s,c_th,SNR_t)

%% Inputs

% L - Number of secondary data
% s - Steering vector
% R_inv - Inverse of the covariance matrix
% xt - Steering vector of the target


%% Output

%p - Approximate value of the PDF of the loss factor evaluated at x

%% Initialization

N = length(s);
M = L-N+1;
x = linspace(0,1,10001);

s_th = 1-c_th;
c = SNR_t*s_th;

p = zeros(1,length(x));

for m = 0:M+1
    c_ = exp(gammaln(N+M-1) - gammaln(N+M-1+m));
    p = p + nchoosek2(M+1,m)*c_*c^m*betapdf_(x,M+1,N+m-1); 
end

p = p.*exp(-c*x);

end

function y = betapdf_(x,n,m)

c = exp(gammaln(n+m) - gammaln(n) - gammaln(m));
y = c*x.^(n-1).*(1-x).^(m-1);

end

function w = nchoosek2(N,K)
    w = exp(gammaln(N+1) - gammaln(K+1) - gammaln(N-K+1));
end




