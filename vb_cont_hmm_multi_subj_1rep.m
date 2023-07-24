function model = vb_cont_hmm_multi_subj_1rep(data,K,maxIter,tol,task,init_option,ss)

%Variational Bayesian Continuous Hidden Markov Model
%
%net = vbhmm(data,K,maxIter,tol);
%
% data (cell with N items, each M*T1 size)
% K - number of states
% maxIter - maximum number of iterations of VB EM (100)
% tol - termination tolerance, prop. change in per-datum likelihood (0.0001)

%best_model is a structure consisting of:
%
% best_model.Wa - state transition Dirichlet counts
% best_model.b - observation emission Wishhart Dist
% best_model.Wpi - initial state prior Dirichlet counts
% best_model.F - F learning curve
%
% Iterates until a proportional change < tol in the per-sequence log
% likelihood, or its iterations of VB.
%
% Srikanth Ryali 10/17/2014
M = size(data{1},1); % no of regions
N = length(data);   % no of subjects
T = zeros(1,N); %length of each dataset
for n = 1:N
    T(n) = size(data{n},2);
end
% Initialise the hyperparameters
alphaa = 1;   % for transition
alphapi = 1;  %for inital prob pi
mi = zeros(M,1);    %prior mean
lambdai = 10^3; %precisionm - large value to force emission posterior means to zero; small values - non-informative
ai = M + 1;      %for Wsihart Distribution W(a,b) - Initial Degrees of Freedom - PAMI-2006 notation
bi = 10^-3 *eye(M); %for Wishart Distribution W(a,b), b = prior precision - PAMI - 2006 notation
% Initialise the pseudo-counts
ua = ones(1,K)*(alphaa/K);
upi = ones(1,K)*(alphapi/K);

s = RandStream.create('mt19937ar', 'Seed', ss);
%RandStream.setDefaultStream(s);
RandStream.setGlobalStream(s);
net = [];
iter = 1;
converged = 0;
flag_violation = 0;
F = [];
while (iter <= maxIter && ~converged && ~flag_violation)
    if(iter == 1 ) % Initialize by Kmeans
        [mp,lambdap,ap,bp,Wa,Wpi] = init_kmeans(data,K,lambdai);
        a = repmat(ai,[K,1]);
        b = repmat(bi,[1,1,K]);
    else
        % M Step
        [mp,lambdap,Wa,Wpi,ap,bp] = vbhmm_mstep(data,phthtpgV1T,phtgV1T,K,ua,upi,ai,bi,mi,lambdai);
    end
    astar = exp(  digamma(Wa) - repmat( digamma(sum(Wa,2)) ,[1 K])  );
    pistar = exp(  digamma(Wpi) - digamma(sum(Wpi,2))  );
    % E Step
    [ind_loglik,phtgV1T,phthtpgV1T] = vbhmm_estep_multi_subj(data,astar',pistar,ap,bp,mp,lambdap); %astar is transposed to match Barber's implementation
    loglik = sum(ind_loglik);
    % Compute F, straight after E Step.
    Fa=0;  
    for kk = 1:K,
        Fa = Fa - kldirichlet(Wa(kk,:),ua);
        %Fb = Fb - klnormal_wishart1(mp,lambdap,ap,bp,mi,lambdai,ai,bi);
    end;
    Fb = -klnormal_wishart2(mp,lambdap,ap,bp,mi,lambdai,ai,bi);
    Fpi = - kldirichlet(Wpi,upi);
    F(iter) = loglik + Fpi + Fa + Fb;
    %Fold = F(iter);
    if (iter >= 2)
        delta_loglik = abs(F(iter) -F(iter-1));
        avg_loglik = (abs(F(iter)) + abs(F(iter-1)) + eps)/2;
        delta_change = delta_loglik / avg_loglik;
        if (F(iter) < F(iter-1))
            converged = 0;
        elseif(F(iter) < F(iter-1) && iter > 1 )
            flag_violation = 1;
            converged = 1;
        elseif (delta_change < tol  && iter > 1)
                 converged = 1;
        end;
    end;
    iter = iter + 1;
end

if converged
    fprintf('seed no = %d \t iterations = %d  \t Converged \t delta_change = %f \t maxF = %f \n',ss, length(F), delta_change, max(F))
elseif flag_violation
    fprintf('seed no = %d \t iterations = %d  \t Violation  \t delta_change = %f \t  maxF = %f \n',ss, length(F), delta_change, max(F))
elseif iter-1 == maxIter
    fprintf('seed no = %d \t iterations = %d  \t Maximum Iters \t delta_change = %f  \t maxF = %f \n',ss, length(F), delta_change, max(F))
end

net.Wa = Wa;
net.mp = mp;
net.lambdap = lambdap;
net.ap = ap;
net.bp = bp;
net.Wpi = Wpi;
net.state_prob_smooth = phtgV1T;
net.F = F;
net.flag_violation = flag_violation;
model = net;
