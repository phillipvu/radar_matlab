function [pd, pf] = off_grid(S,S_est,ep)

if isempty(S_est) == 1
    pf = 0;
    pd = 0;
    return
end

for i = 1:size(S,1)
    
    [T,I] = min(sum(abs(S_est - repmat(S(i,:),[size(S_est,1),1])).^2,2));
    if sqrt(T) < ep
        S_est(I) = inf;
    end
        
end

pd = sum(isinf(S_est(:,1)));
pf = length(S_est(:,1))-pd;

end