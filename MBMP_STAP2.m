function S_ = MBMP_STAP2(Y,A,D,e,Z)

%% Inputs

% Y - Measurement vector
% A - Dictionary matrix
% D - Branch vector
% e - threshold parameter
% Z - Secondary data

%% Outputs

% S_ - Set of targets

%% Initialization

[N,L] = size(Z);

% R = Z*Z'/L;
R = sqrtm(Z);
Y = R*Y;
A = R*A;

D = [D ones(1,size(A,2)-length(D))];

S_ = [];
I_t = 1;
L = 1;

for s = 1:length(e)
   S_(s).S = [];
   S_(s).e = e(s);
end

e = 1-e;

%% Branching

B = [];
B.S = [];
B.proj = eye(size(Y,1));
B.fro = 2*norm(Y,'fro')^2;
B.lvl = 0;

B_ = B.fro;

x = chi2inv(e,2*N);

for f = I_t:length(x)

    if B_ <= x(f)
        I_t = I_t+1;
    else
        break
    end
end

if I_t > length(x)
    return
end

for d = 1:length(D)

for i = 1:L

A_ = B(i).proj*A;
normA = sum(abs(A_).^2,1)';
normA(normA < 1e-12) = inf;

[~,I] = sort(abs(A_'*Y).^2./normA,'descend');

B_temp = tag(Y,A,B(i),I(1:D(d)));
B(end+1:end+length(B_temp)) = B_temp; 

end

B(1:L) = [];
S_k = [];

for l = 1:length(B)
    S_k(l,:) = sort(B(l).S);
end
   
[~,U,~] = unique(S_k,'rows');
B = B(sort(U));

[~,I] = min([B.fro]);
B_ = B(I).fro;

x = chi2inv(e,2*(N-B(I).lvl));

for f = I_t:length(x)

    if B_ <= x(f)
        S_(f).S = B(I).S;
        I_t = I_t+1;
    else
        break
    end
end

if I_t > length(x)
    return
end

L = length(B);

end

end

    
function B = tag(Y,A,B_parent,I)
B = [];  
for i = 1:length(I)
    
    B(end+1).S = [I(i) B_parent.S];
    proj_ = B_parent.proj*A(:,I(i));
    B(end).proj = B_parent.proj - proj_*proj_'/norm(proj_)^2; 
    B(end).fro = 2*norm(B(end).proj*Y,'fro')^2;
    B(end).lvl = length(B(end).S);    
        
end
    
end
