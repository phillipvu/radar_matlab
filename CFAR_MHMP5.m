function S_ = CFAR_MHMP5(Y,zr,zt,P,g_,n,D,t,Z)

%% Inputs

% Y - Received data
% zr - Sequence of receieve elements
% zt - Sequence of transmit elements
% P - Number of pulses used
% g_ - set of grid points used in each dimension
% n - Desired granularity
% D - Branch vector
% t - Desired false alarm probability
% Z - Set of secondary data

%% Output

% S_ - Set of angle-Doppler points that correspond to targets

I_t = 1;
[N,L] = size(Z);
del_ = g_(2)-g_(1);
[A,AD] = generate_dictionary2(zr,zt,g_,g_,P);
S_ = [];

for i = 1:length(t)
    S_(i).AD = [];
end

R = Z*Z';
R = sqrtm(inv(R));
Y = R*Y;
A = R*A;

D = [D ones(1,size(A,2)-length(D))];
x = 1-(1-t).^(1/size(A,2));
g = (x.^(-1/(L-N+2))-1)*(L+1)/(L-N+1);

for k = 1:length(g_)^2
    del = del_;
    A_ = A;
    AD_ = AD;
 
    B = MBMP(Y,A_,D(1:k));
    [~,I] = min([B.fro]);
    S = unique([B.S]);

    A_2 = [];
    AD_2 = [];

    for s = 1:length(S)

        gp = AD_(S(s),1)-del:n:AD_(S(s),1)+del;
        gp(gp < -1) = [];
        gp(gp > 1) = [];

        ga = AD_(S(s),2)-del:n:AD_(S(s),2)+del;
        ga(ga < -1) = [];
        ga(ga > 1) = [];


        for p = 1:length(gp)
            for a = 1:length(ga)
                AD_2 = [AD_2;gp(p) ga(a)];
            end
        end

    end

    [~,IA,~] = unique(AD_2,'rows');  
    AD_ = AD_2(IA,:);
    A_2 = generate_dictionary3(zr,zt,AD_2,P);
    A_ = R*A_2(:,IA);
            
    B = MBMP(Y,A_,D(1:k));
    [~,I] = min([B.fro]);
    S = B(I).S;
    
    T = ASD(Y,A_(:,S));
        
    for f = I_t:length(g)

        if T > g(f)
            S_(f).AD = AD_(S,:);
        else
            I_t = f+1; 
        end

    end

    if I_t > length(g)
        return
    end
        
end
    


end



function B = MBMP(Y,A,D) 

%% Inputs

% Y - Received data
% A - Set of steering vectors
% D - Branch vector

%% Output

% B - Nodes at the end of the tree

%% Initialization

S_ = [];

%% Branching

B = [];
B.S = [];
B.proj = eye(size(Y,1));
B.fro = norm(Y,'fro')^2;
L = 1;


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
L = length(B);

end

end
    
function B = tag(Y,A,B_parent,I)
B = [];  
for i = 1:length(I)
    
    B(end+1).S = [I(i) B_parent.S];
    proj_ = B_parent.proj*A(:,I(i));
    B(end).proj = B_parent.proj - proj_*proj_'/norm(proj_)^2; 
    B(end).fro = norm(B(end).proj*Y,'fro')^2;
        
end
    
end




function T = ASD(Y,A)
    
    T = sum(abs(A\Y).^2,2);
    A_ = diag(inv(A'*A));
    
    T = abs(T./A_);


end



