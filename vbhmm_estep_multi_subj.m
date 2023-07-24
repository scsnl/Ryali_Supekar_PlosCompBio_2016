function [loglik,phtgV1T,phthtpgV1T] = vbhmm_estep_multi_subj(v,stran,sprior,ap,bp,mp,lambdap)

%EM training of a Switching AR model
% Inputs:
% v :  timeseries Regions * Time (M*T)
% stran : transition probabilities
%sprior - initial probablity
%a,b - Wishart parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SFC.a = ap;
SFC.b = bp;
SFC.mp = mp;
SFC.lambdap = lambdap;
parfor n = 1:length(v)
    % Estep
    [logalpha,loglik1]=HMMforwardSFC(v{n},stran,sprior,SFC);
    logbeta=HMMbackwardSFC(v{n},stran,SFC);
    [phtgV1T1,phthtpgV1T1]=HMMsmoothSFC(logalpha,logbeta,SFC,stran,v{n});
    
    % Transpose joint probs phthtpgV1T to match the notation of Murphy, Beal & PAMI-2006
    for t = 1:size(phthtpgV1T1,3)
        phthtpgV1T1(:,:,t) = phthtpgV1T1(:,:,t)';
    end
    loglik(n) = loglik1;
    phtgV1T{n} = phtgV1T1;
    phthtpgV1T{n} = phthtpgV1T1;
end










