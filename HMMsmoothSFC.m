function [phtgV1T,phthtpgV1T1]=HMMsmoothSFC(logalpha,logbeta,SFC,phghm,v)
%HMMSMOOTHSFC Switching Autoregressive HMM smoothing
% [phtgV1T,phthtpgV1T]=HMMsmoothSFC(logalpha,logbeta,a,sigma2,phghm,v,Tskip)
% return the smoothed pointwise posterior p(h(t)|v(1:T)) and pairwise smoothed posterior p(h(t),h(t+1)|v(1:T)).
%
% Inputs:
% logalpha : log alpha messages (see HMMforwardSFC.m)
% logbeta : log beta messages (see HMMbackwardSFC.m)
% % phghm : state (switch) transition matrix p(h(t)|h(t-1))
% v : observations
% Tskip : the number of timesteps to skip before a switch update is allowed
%
% Outputs:
% phtgV1T : smoothed posterior p(h(t)|v(1:T))
% phthtpgV1T  : smoothed posterior p(h(t),h(t+1)|v(1:T))
% See also HMMforwardSFC.m, HMMbackwardSFC.m, demoSFClearn.m
T = size(v,2);  %length of time series
H = length(SFC.a);  % # of states
M = size(v,1); %# of regions
% smoothed posteriors: pointwise marginals:
for t=1:T
    logphtgV1T(:,t)=logalpha(:,t)+logbeta(:,t); % alpha-beta approach
    phtgV1T(:,t)=condexp(logphtgV1T(:,t));
end
% smoothed posteriors: pairwise marginals p(h(t),h(t+1)|v(1:T)):
for t=2:T
    atmp=condexp(logalpha(:,t-1));
    btmp=condexp(logbeta(:,t));
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
        phatvgh1(h) = exp(term1 + term2 + term3 + term4 + term5) + eps;
    end
    phatvgh1=condexp(phatvgh1);
    phghmt=phghm;
    ctmp1 = repmat(atmp,1,H).*phghmt'.*repmat(phatvgh1'.*btmp',H,1)+eps; % two timestep potential
    phthtpgV1T1(:,:,t-1)=ctmp1./sum(sum(ctmp1));
end