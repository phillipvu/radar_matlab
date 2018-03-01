function T = AMF_CFAR(y,s,I,Z)

%% Inputs

% y - Measurement vector
% s - Signal vector
% I - Interference subspace
% Z - Secondary data

%% Output

%T - Output of the CFAR detector 

R_inv = inv(Z*Z');

if isempty(I) == 1
    d = s;
else
    N = length(y);
    d = (eye(N) - I*inv(I'*R_inv*I)*I'*R_inv)*s;
end

R_inv = inv(Z*Z');
T = abs(d'*R_inv*y)^2/abs(d'*R_inv*d);

end