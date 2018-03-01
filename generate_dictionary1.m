function A = generate_dictionary(z,t,g,p)

%% Inputs

% z - Sequence of receive element positions, assumed to be a column vector
% t - Sequence of transmitter element positions, assumed to be a column vector set to '0' if not MIMO array
% g - number of grid points to test
% p - number of pulses set to 0 for angle only dimension 

%% Outputs

% A - Dictionary matrix, columns are steering vectors on the angle-Doppler map (or angle dimension only)
Fc = 1e9; % carrier frequency
Fs = 40e6; % sampling frequency
phase_noise_freq = [ 1e3, 10e3, 100e3, 1e6, 10e6 ]; % Offset From Carrier
phase_noise_power = [ -84, -100, -96, -109, -122 ]; % Phase Noise power

A = [];
N = length(z)*length(t);

if p > 0

L = 0.5*(0:p-1)';

for i = 1:length(g)
    b = exp(1i*2*pi*L*g(i)*Fc/Fs);
     b1 = add_phase_noise( b, Fs, phase_noise_freq, phase_noise_power );
    for j = 1:length(g)
        at = exp(1i*2*pi*t*g(j)*Fc/Fs);
         at1 = add_phase_noise( at, Fs, phase_noise_freq, phase_noise_power );
        ar = exp(1i*2*pi*z*g(j)*Fc/Fs);
         ar1 = add_phase_noise( ar, Fs, phase_noise_freq, phase_noise_power );
        c = kron(at1,ar1);
      
        A(:,end+1) = kron(b1,c)/sqrt(N*p);
    end
end

else
    
    for j = 1:length(g)
%         at = exp(1i*2*pi*t*g(j));
%        % at1 = add_phase_noise( at, Fs, phase_noise_freq, phase_noise_power );
%         ar = exp(1i*2*pi*z*g(j));
       % ar1 = add_phase_noise( ar, Fs, phase_noise_freq, phase_noise_power );
         at = exp(1i*2*pi*t*g(j)*Fc/Fs);
         at1 = add_phase_noise( at, Fs, phase_noise_freq, phase_noise_power );
        ar = exp(1i*2*pi*z*g(j)*Fc/Fs);
         ar1 = add_phase_noise( ar, Fs, phase_noise_freq, phase_noise_power );
         
        A(:,end+1) = kron(at1,ar1)/sqrt(N);
    end
        
end


end