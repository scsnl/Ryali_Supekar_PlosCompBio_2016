function  [ap,bp] = update_Wishartparams(data,ai,bi,phtgV1T)


% Update the parameters for the Wishart Distribution

w = sum(phtgV1T,2);
K = size(phtgV1T,1);   % # of states
M = size(data,1);      % # of regions
S = zeros(M,M,K);      % sufficient stats for Ri - precision
for t = 1:size(data,2)
    for k = 1:K
        S(:,:,k) = S(:,:,k) + phtgV1T(k,t) * data(:,t)*data(:,t)';
    end
end
 ao = repmat(ai,[K,1]);
 bo = repmat(bi,[1,1,K]);
 ap = ao + w;
 bp = bo + S;


