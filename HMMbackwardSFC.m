function logbeta1 = HMMbackwardSFC(v,phghm,SFC)
%HMMBACKWARDSFC Backward Pass (beta method) for the Switching Autoregressive HMM
% logbeta=HMMbackwardSFC(v,phghm,a,sigma2,Tskip)
%
% Inputs:
% v : observations
% phghm : state (switch) transition matrix p(h(t)|h(t-1))
% Tskip : the number of timesteps to skip before a switch update is allowed
%
% Outputs:
% logbeta: log backward messages log p(v(t+1:T)|h(t),v(t-L+1:t))
% See also HMMforwardSFC.m and demoSFClearn.m
T = size(v,2);  %length of time series
H = length(SFC.a);  % # of states
M = size(v,1); %# of regions
% logbeta recursion
logbeta1 = zeros(H,T);
logbeta1(:,T) = zeros(H,1);
for t=T:-1:2
    phatvgh1 = zeros(H,1);
   for h = 1:H
        %for Wishart distribution
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
        phatvgh1(h) = exp(term1 + term2 + term3 + term4 + term5)+eps; % Eq 44
   end
    logbeta1(:,t-1)=logsumexp(repmat(logbeta1(:,t),1,H),repmat(phatvgh1,1,H).*phghm);
end