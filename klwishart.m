function res =  klwishart(ap,bp,ai,bi)

%KL(q||p), q ~ Wishart(ap,bp), p ~ Wishart(ai,bi)

res = 0;
M = size(bp,1);
K = size(bp,3);
for k = 1:K
    %compute log_Zq
    gamma_args = repmat(ap(k) + 1,1,M)-(1:M) ;
    log_Zq = -0.5*ap(k)*log(det(0.5*bp(:,:,k))); +  sum(log(gamma(0.5*gamma_args)));
    %Compute log_Zp
    gamma_args = repmat(ai(k) + 1,1,M)-(1:M) ;
    log_Zp = -0.5*ai(k)*log(det(0.5*bi(:,:,k))); +  sum(log(gamma(0.5*gamma_args)));
    log_ratio = log_Zp-log_Zq;
    %Compute E(log(|R|))
    digamma_args = repmat(ap(k) + 1,1,M)-(1:M) ;
    exp_log_detR_q = -log(det(0.5*bp(:,:,k))) + sum(digamma(0.5*digamma_args));
    digamma_args = repmat(ai(k) + 1,1,M)-(1:M) ;
    exp_log_detR_p = -log(det(0.5*bi(:,:,k))) + sum(digamma(0.5*digamma_args));
    %res = res + 0.5*(ai(k) - ap(k))* exp_log_detR - 0.5*ap(k)*M + 0.5*ap(k)*trace(bi(:,:,k)*pinv(bp(:,:,k))) + log_ratio;
    res = res + 0.5*(ap(k)-M-1)*exp_log_detR_q -0.5*(ai(k)-M-1)*exp_log_detR_p - 0.5*ap(k)*M + 0.5*ap(k)*trace(bi(:,:,k)*pinv(bp(:,:,k))) + log_ratio;
end