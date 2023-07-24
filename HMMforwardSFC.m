function [logalpha1,loglik1] = HMMforwardSFC(v,phghm,ph1,SFC)
%HMMFORWARDSFC Switching Autoregressive HMM with switches updated only every Tskip timesteps
% [logalpha,loglik]=HMMforwardSFC(v,phghm,ph1,a,sigma2,Tskip)
%
% Inputs:
% v : observations
% phghm : state (switch) transition matrix p(h(t)|h(t-1))
% ph1 : prior state distribution
% % Tskip : the number of timesteps to skip before a switch update is allowed
%
% Outputs:
% logalpha : log forward messages
% loglik : sequence log likelihood log p(v(1:T))
% See also HMMbackwardSFC.m and demoSFClearn.m
T = size(v,2);  %length of time series
H = length(ph1);  % # of states
M = size(v,1); %# of regions
logalpha1 = zeros(H,T);
% logalpha recursion:
for h = 1:H
    %for Normal - Wishart distribution
    b = SFC.b(:,:,h);
    a = SFC.a(h);
    lambdap = SFC.lambdap(h);
    mp = SFC.mp(:,h);
    term1 = -M/2 * log(2*pi);
    term2 = -0.5*log(det(0.5*b));
    digamma_args = repmat(a + 1,1,M)-(1:M) ;
    term3 = 0.5*sum(digamma(0.5*digamma_args));
    term4 = -0.5 *a* (v(:,1)-mp)'*pinv(b)*(v(:,1)-mp);
    term5 = -0.5*M/lambdap;
    logalpha1(h,1) = term1 + term2 + term3 + term4 + term5 + log(ph1(h)); % Eq 44
end
for t = 2:T
    phatvgh1 = zeros(H,1);
    for h = 1:H
        %for Normal- Wishart distribution
        b = SFC.b(:,:,h);
        a = SFC.a(h);
        lambdap = SFC.lambdap(h);
        mp = SFC.mp(:,h);
        term1 = -M/2 * log(2*pi);
        term2 = -0.5*log(det(0.5*b));
        digamma_args = repmat(a + 1,1,M)-(1:M) ;
        term3 = 0.5*sum(digamma(0.5*digamma_args));
        term4 = -0.5*a*(v(:,t)-mp)'*pinv(b)*(v(:,t)-mp);
        term5 = -0.5*M/lambdap;
        phatvgh1(h) = exp(term1 + term2 + term3 + term4 + term5) + eps; % Eq 44
    end
    logalpha1(:,t)=logsumexp(repmat(logalpha1(:,t-1),1,H),repmat(phatvgh1',H,1).*phghm');
end
loglik1 = logsumexp(logalpha1(:,T),ones(H,1)); % log likelihood