function best_model = vb_cont_hmm(data,K,maxIter,tol, no_repetitions,task,init_option)

%Variational Bayesian Continuous Hidden Markov Model
%
%net = vbhmm(data,K,maxIter,tol);
%
% data (M*T)
% K - number of states
% maxIter - maximum number of iterations of VB EM (100)
% tol - termination tolerance, prop. change in per-datum likelihood (0.0001)

%net is a structure consisting of:
%
% net.Wa - state transition Dirichlet counts
% net.b - observation emission Wishhart Dist
% net.Wpi - initial state prior Dirichlet counts
% net.F - F learning curve
%
% Iterates until a proportional change < tol in the per-sequence log
% likelihood, or its iterations of VB.
%
% Srikanth Ryali 10/17/2014


[M,total_length] = size(data); % no of regions, length of data
% Initialise the hyperparameters
alphaa = 1;   % for transition
alphapi = 1;  %for inital prob pi
ai = M;      %for Wsihart Distribution W(a,b) - Initial Degrees of Freedom - PAMI-2006 notation
bi = 10^-3 *eye(M); %for Wishart Distribution W(a,b), b = prior precision - PAMI - 2006 notation
% Initialise the pseudo-counts
ua = ones(1,K)*(alphaa/K);
upi = ones(1,K)*(alphapi/K);
if init_option == 2
    no_repetitions = 1;
end
for repetition = 1: no_repetitions
    iter = 1;
    converged = 0;
    flag_converged = 0; flag_violation = 0;
    F = [];
    while (iter <= maxIter && ~converged)
        if (iter==1 && init_option == 1)
            % Initialization
            % Pick an HMM from the prior to initialize the counts
            wa = [];
            for k=1:K, % loop over hidden states
                wa(k,:) = dirrnd(ua,1)*total_length;
            end;
            wpi = dirrnd(upi,1);
            Wa = wa + repmat(ua,[K 1]);
            Wpi = wpi + upi;
            ap = repmat(ai,[K,1]);
            bp = repmat(bi,[1,1,K]);
            a = ap; b = bp;
        elseif(iter == 1 && init_option == 2)
            [ap,bp,Wa,Wpi] = init_vbhmm(data,task);
            a = repmat(ai,[K,1]);
            b = repmat(bi,[1,1,K]);
        elseif(iter == 1 && init_option == 3)
            [ap,bp,Wa,Wpi] = init_kmeans(data,K);
            a = repmat(ai,[K,1]);
            b = repmat(bi,[1,1,K]);
        else
            % M Step
            %compute parameters for Dirichlet
            wa = sum(phthtpgV1T,3); % q(s_t = i, s_t+1 = j)
            Wa = wa + repmat(ua,[K 1]); % posterior is data counts plus prior.
            wpi = phtgV1T(:,1)';     % q(s1 = i)
            Wpi = wpi + upi; % posterior is data counts plus prior
            %update Wishart parameters
            [ap,bp] = update_Wishartparams(data,ai,bi,phtgV1T);
        end
        astar = exp(  digamma(Wa) - repmat( digamma(sum(Wa,2)) ,[1 K])  );
        pistar = exp(  digamma(Wpi) - digamma(sum(Wpi,2))  );
        % E Step
        [loglik,phtgV1T,phthtpgV1T] = vbhmm_estep(data,astar',pistar,ap,bp);
        % Compute F, straight after E Step.
        Fa=0; Fb=0; Fpi=0;
        for kk = 1:K,
            Fa = Fa - kldirichlet(Wa(kk,:),ua);
            Fb = Fb - klwishart(ap,bp,a,b);
        end;
        Fpi = Fpi - kldirichlet(Wpi,upi);
        F(iter) = Fa + Fb + Fpi + loglik;
        %Fold = F(iter);
        if (iter >= 2)
            delta_loglik = abs(F(iter) -F(iter-1));
            avg_loglik = (abs(F(iter)) + abs(F(iter-1)) + eps)/2;
            delta_change = delta_loglik / avg_loglik;
            if (F(iter) < F(iter-1) && iter <= 10)
                converged = 0;
            elseif(F(iter) < F(iter-1) && iter > 10 )
                flag_violation = 1;
                converged = 1;
            elseif (delta_change < tol  && iter >= 10)
                flag_converged = 1;
                converged = 1;
            end;
        end;
        iter = iter + 1;
    end;
    if flag_converged
        fprintf('repetition no = %d \t iterations = %d  \t Converged \t delta_change = %f \t maxF = %f \n',repetition, length(F),delta_change,max(F))
     elseif flag_violation
        fprintf('repetition no = %d \t iterations = %d  \t Violation  \t delta_change = %f \t  maxF = %f \n',repetition, length(F),delta_change,max(F))
     elseif iter-1 == maxIter
        fprintf('repetition no = %d \t iterations = %d  \t Maximum Iters \t delta_change = %f  \t maxF = %f \n',repetition, length(F),delta_change,Max(F))
    end
    net.Wa = Wa;
    net.ap = ap;
    net.bp = bp;
    net.Wpi = Wpi;
    net.state_prob_smooth = phtgV1T;
    net.F = F;
    if repetition == 1
        best_model = net;
    elseif (max(best_model.F) < max(F)) 
        best_model = net;
    end
end

    
