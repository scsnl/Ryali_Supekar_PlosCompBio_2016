function [mp,lambdap,Wa,Wpi,ap,bp] = vbhmm_mstep(data,phthtpgV1T,phtgV1T,K,ua,upi,ai,bi,mi,lambdai)

wa = zeros(K,K); wpi = zeros(1,K);
for n = 1:length(data)
    phthtpgV1T1 = phthtpgV1T{n}; phtgV1T1 = phtgV1T{n};
    wa = wa + sum(phthtpgV1T1,3); % q(s_t = i, s_t+1 = j)
    wpi = wpi + phtgV1T1(:,1)'; % q(s1 = i)
end
Wa = wa + repmat(ua,[K 1]); % posterior is data counts plus prior.
Wpi = wpi + upi; % posterior is data counts plus prior

% Update the parameters for the Normal-Wishart Distribution
w = zeros(K,1);
M = size(data{1},1);  % # of regions
S = zeros(M,M,K);      % sufficient stats for Ri - precision
xb = zeros(M,K);        %sufficient stats for posterior mean
for n = 1:length(data)
    phtgV1T1 = phtgV1T{n};
    w = w + sum(phtgV1T1,2); % Eq 32
end
for k = 1:K
    for n = 1:length(data)
        phtgV1T1 = phtgV1T{n};
        data1 = data{n};
        wk = ones(M,1) * phtgV1T1(k,:);
        xb(:,k) = xb(:,k) + sum((wk.* data1),2); %Eq 33
    end
    xb(:,k) = xb(:,k)./w(k); % Eq 33 : normalize by w(k)
end
for n = 1:length(data)
    phtgV1T1 = phtgV1T{n};
    data1 = data{n};
    for t = 1:size(data1,2)
        for k = 1:K
            S(:,:,k) = S(:,:,k) + phtgV1T1(k,t) * (data1(:,t)-xb(:,k))*(data1(:,t)-xb(:,k))'; %Eq 34
        end
    end
end
ao = repmat(ai,[K,1]);
ap = ao + w; %Eq 35
lambdap = ones(K,1)*lambdai + w; %Eq 37
mp = zeros(M,K); bp = zeros(M,M,K);
for k = 1:K
    bp(:,:,k) = bi + S(:,:,k) + (lambdai * w(k) * (mi - xb(:,k)) * (mi - xb(:,k))')/lambdap(k); %Eq 36
    mp(:,k) = (lambdai * mi + w(k)*xb(:,k))./lambdap(k); % Eq 38
end





