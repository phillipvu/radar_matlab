function A = generate_dictionary3(z,t,AD,p)

%% Inputs

% z - Sequence of receive element positions, assumed to be a column vector
% t - Sequence of transmitter element positions, assumed to be a column vector set to '0' if not MIMO array
% AD - Set of angle-Doppler points to test
% p - number of pulses set to 0 for angle only dimension 

%% Outputs

% A - Dictionary matrix, columns are steering vectors on the angle-Doppler map (or angle dimension only)
% S - Angle-Doppler pairs that correspond to the dictionary matrix A

A = [];
N = length(z)*length(t);
L = 0.5*(0:p-1)';

for i = 1:size(AD,1)
    b = exp(1i*2*pi*L*AD(i,1));

    at = exp(1i*2*pi*t*AD(i,2));
    ar = exp(1i*2*pi*z*AD(i,2));
    c = kron(at,ar);

    A(:,end+1) = kron(b,c)/sqrt(N*p);

end




end