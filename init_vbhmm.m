function [ap,bp,Counter,sprior] = init_vbhmm(v,states)

S = length(unique(states)); % # of states
[M,T] = size(v);         % # of Regions , # of time samples
phtgV1T = zeros(S,T);
sprior = zeros(1,S);
for s = 1:S
    ix = states == s;
    sprior(s) = sum(ix);
end
% Estimate Covariance for each state
ap = zeros(S,1); bp = zeros(M,M,S);
for s = 1:S
    ix = find(states == s);
    ap(s) = length(ix);
    y = v(:,ix)';
    bp(:,:,s) = cov(y) + eps.*eye(M);
end
%Estimate Transition Probability Matrix
Counter = zeros(S);
for t1 = 1:T-1
    t2 = t1+1;
    m = states(t2); n = states(t1);
    Counter(m,n) = Counter(m,n) + 1;
end