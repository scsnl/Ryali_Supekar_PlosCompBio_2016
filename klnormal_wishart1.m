function res =  klnormal_wishart1(mp,lambdap,ap,bp,mi,lambdai,ai,bi)

%KL(q||p), q ~ Normal-Wishart(mp,lambdap,ap,bp), p ~ Normal-Wishart(mi,lambdai,ai,bi)

res = 0;
M = size(bp,1);
K = size(bp,3);
for k = 1:K
    %compute log_Zq
    gamma_args = repmat(ap(k) + 1,1,M)-(1:M) ;
    log_Zq = 0.5*ap(k)*M*log(2)-0.5*ap(k)*log(det(bp(:,:,k))); +  sum(gammaln(0.5*gamma_args));
    %Compute log_Zp
    gamma_args = repmat(ai + 1,1,M)-(1:M) ;
    log_Zp = 0.5*ai*M*log(2)-0.5*ai*log(det(bi)); +  sum(gammaln(0.5*gamma_args));
    log_ratio = log_Zp-log_Zq;
    %Compute E(log(|R|))
    digamma_args = repmat(ap(k) + 1,1,M)-(1:M) ;
    exp_log_detR_q = -log(det(bp(:,:,k))) + sum(digamma(0.5*digamma_args));
    digamma_args = repmat(ai + 1,1,M)-(1:M) ;
    exp_log_detR_p = -log(det(bi)) + sum(digamma(0.5*digamma_args));
    %res = res + 0.5*(ai(k) - ap(k))* exp_log_detR - 0.5*ap(k)*M + 0.5*ap(k)*trace(bi(:,:,k)*pinv(bp(:,:,k))) + log_ratio;
    %KL - Wishart
    res = res + 0.5*(ap(k)-M-1)*exp_log_detR_q -0.5*(ai-M-1)*exp_log_detR_p - 0.5*ap(k)*M + 0.5*ap(k)*trace(bi*pinv(bp(:,:,k))) + log_ratio;
    %Additional Terms
    res = res + 0.5*(M*log(lambdap(k)/lambdai) + M*lambdap(k)/lambdai - M + lambdai * (mp(:,k)-mi)'*ap(k)*pinv(bp(:,:,k)) * (mp(:,k)-mi));   
end