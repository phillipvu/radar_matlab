function C = generate_clutter(z,t,n,p)

%% Inputs

% z - Sequence of receive element positions, assumed to be a column vector
% t - Sequence of transmitter element positions, assumed to be a column vector set to '0' if not MIMO array
% p - number of pulses set to 

%% Outputs

% C - Clutter matrix, columns are steering vectors on the angle-Doppler map (or angle dimension only)

C = [];
N = length(z)*length(t);
g = linspace(-1,1,n);
L = 0.5*(0:p-1)';

for i = 1:length(g)
    b = exp(1i*2*pi*L*g(i));

        at = exp(1i*2*pi*t*g(i));
        ar = exp(1i*2*pi*z*g(i));
        
        C(:,end+1) = kron(b,kron(at,ar))/sqrt(N*p);
end


end