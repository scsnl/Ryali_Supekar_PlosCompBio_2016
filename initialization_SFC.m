function [SFC,stran,sprior] = initialization_SFC(v,states)

S = length(unique(states)); % # of states
[M,T] = size(v);         % # of Regions , # of time samples
phtgV1T = zeros(S,T);
sprior = zeros(S,1);
for s = 1:S
    ix = states == s;
    sprior(s) = sum(ix)/T;
end
% Estimate Covariance for each state
for s = 1:S
    ix = states == s;
    y = v(:,ix)';
    %SFC(s).Q = cov(y);
    SFC(s).Q = cov(y) + 10^-10.*eye(M);
    SFC(s).mean = mean(y,1)';
end
%Estimate Transition Probability Matrix
Counter = zeros(S);
for t1 = 1:T-1
    t2 = t1+1;
    m = states(t2); n = states(t1);
    Counter(m,n) = Counter(m,n) + 1;
end
stran = condp(Counter);