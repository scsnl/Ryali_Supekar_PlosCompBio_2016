function [loglik,phtgV1T,phthtpgV1T] = vbhmm_estep(v,stran,sprior,a,b)

%EM training of a Switching AR model
% Inputs:
% v :  timeseries Regions * Time (M*T)
% stran : transition probabilities
%sprior - initial probablity
%a,b - Wishart parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SFC.a = a;
SFC.b = b;

% Estep
[logalpha,loglik]=HMMforwardSFC(v,stran,sprior,SFC);
logbeta=HMMbackwardSFC(v,stran,SFC);
[phtgV1T,phthtpgV1T]=HMMsmoothSFC(logalpha,logbeta,SFC,stran,v);

% Transpose joint probs phthtpgV1T to match the notation of Murphy, Beal & PAMI-2006
for t = 1:size(phthtpgV1T,3)
    phthtpgV1T(:,:,t) = phthtpgV1T(:,:,t)';
end







