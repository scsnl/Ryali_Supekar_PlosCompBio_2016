function [mp,lambdap,ap,bp,Counter,sprior] = init_rand(v,K,lambdai)

data = [];
data_un = [];
%%%%% Standardize data %%%%%%%%%%%%%%
for n = 1:length(v)
    x = v{n};
    data_un = [data_un x];
    for t = 1:size(v{n},2)
        x(:,t) = x(:,t)/norm(x(:,t));
    end
    data = [data x];
end
states = randi([1,K],length(data),1)';
%S = length(unique(states)); % # of states
S = K;
sprior = zeros(1,S);
for s = 1:S
    ix = states == s;
    sprior(s) = sum(ix);
end
[M,T] = size(data);
%estimate mean & lambdap for each state
mp = zeros(M,S);
lambdap = zeros(1,S);
for s = 1:S
    ix = find(states == s);
    lambdap(s) = lambdai + length(ix)/T;
    y = data_un(:,ix);
    mp(:,s) = mean(y,2);
end
% Estimate Covariance for each state
ap = zeros(S,1); bp = zeros(M,M,S);
for s = 1:S
    ix = find(states == s);
    ap(s) = length(ix);
    y = data_un(:,ix)';
    bp(:,:,s) = cov(y); + eps.*eye(M);
end
%Estimate Transition Probability Matrix
Counter = zeros(S);
for t1 = 1:T-1
    t2 = t1+1;
    m = states(t2); n = states(t1);
    Counter(m,n) = Counter(m,n) + 1;
end
Counter = Counter';