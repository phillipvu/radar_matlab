function S_ = MBMP_STAP(Y,A,D,t,Z)

%% Inputs

% Y - Measurement vector
% A - Dictionary matrix
% D - Branch vector
% t - threshold parameter
% Z - Secondary data

%% Outputs

% S_ - Set of targets

%% Initialization

[N,L] = size(Z);
R = Z*Z';
R = sqrtm(inv(R));
Y = R*Y;
A = R*A;


x = (1-t).^(1/size(A,2));
g = finv(x,2,2*N);
S_ = [];
I_t = length(t);
L = 1;

for s = 1:length(t)
   S_(s).S = [];
   S_(s).t = t(s);
end

%% Branching

B = [];
B.S = [];
B.proj = eye(size(Y,1));
B.fro = norm(Y,'fro')^2;
B.lvl = 0;
B.T = [];

B_ = B;


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
B_(end+1) = B(I);
L = length(B);

end

for s = length(B_):-1:2

% g = finv(x,2,2*(N-B_(s).lvl));    
I_ = I_t;
    
for f = I_:-1:1
    
    if min([B_(s).T]) >= g(f)
        S_(f).S = [B_(s).S];
        I_t = I_t-1;
    end   
        
end

if I_t == 0
    return
end

end

end

    
function B = tag(Y,A,B_parent,I)
B = [];  
for i = 1:length(I)
    
    B(end+1).S = [I(i) B_parent.S];
    proj_ = B_parent.proj*A(:,I(i));
    B(end).proj = B_parent.proj - proj_*proj_'/norm(proj_)^2; 
    B(end).fro = norm(B(end).proj*Y,'fro')^2;
    B(end).lvl = length(B(end).S);    
    B(end).T = GLRT(Y,A(:,B(end).S),B(end).fro);
        
end
    
end

function T = GLRT(Y,A,fro)
    
    [N,S] = size(A);
        
    T = sum(abs(A\Y).^2,2);
    A_ = diag(inv(A'*A));
    T = (N-S)*abs(T./A_)/fro;

end
