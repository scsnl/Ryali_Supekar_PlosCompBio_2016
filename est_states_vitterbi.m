function [states] =  est_states_vitterbi(data,model)

% Vitterbi Algo


%%     transition distribution p(h(t)|h(t-1)
Wa = model.Wa;
Wa = Wa';
K = size(Wa,1);
Aest = Wa./repmat(sum(Wa,1),K,1); %transition distribution p(h(t)|h(t-1))
Aest = Aest';
%%%%    estimeted pi
Wpi = model.Wpi;
piest = Wpi./sum(Wpi);
% Wishart Parameters
ap = model.ap;
bp = model.bp;
%% Most likely state by Vitterbi Algo
%%%%       emission distribution p(v(t)|h(t)) - Multivariate distribution
states = [];
for n = 1:length(data)
    data1 = data{n};
    [M,T] = size(data1);
    pvgh = zeros(K,T);
    for k = 1:K
        Cov_mat =  bp(:,:,k)/ap(k);
        mu = model.mp(:,k)';
        for t = 1:size(data1,2)
            pvgh(k,t) = mvnpdf(data1(:,t)',mu,Cov_mat);
        end
    end
    states = [states viterbi_path(piest, Aest',pvgh)];
end

