clear all
close all
clc

K = 25; % # of states
maxIter = 100; % maximum number of allowed iterations
tol = 10^-3;   % tolerance for convergence
no_repetitions = 3; % Run vbhmm 100 times and choose the best model
init_option = 3; %1-random, 2-task,3-kmeans
%%%%% simulate data %%%%%%%%
load('stim_design.mat')
std_noise = 0.1;
%load data
%   for n = 1:5
%       [data1,task1] = sim_corr_data(stim_design,std_noise);
%       data{n} = data1; 
%      task{n} = stim_design + 1;
%   end
load simdata_nomean
% VB-HMM
matlabpool local 8
[Model]= vb_cont_hmm_multi_subj_test(data,K,maxIter,tol,no_repetitions,task,init_option);
matlabpool close
for repetition = 1:length(Model)
    F = Model{repetition}.F;
     if repetition == 1
         model = Model{repetition};
     elseif (max(model.F) < max(F))
         model = Model{repetition};
     end
end
%save(fname,'model','task','data')
no_subjs = length(model.state_prob_smooth);
no_obs = size(model.state_prob_smooth{1},2);
[est_states] =  est_states_vitterbi(data,model);
est_states = reshape(est_states',[no_obs,no_subjs]);
state_prob_smooth = [];
for n = 1:no_subjs
    state_prob_smooth = [state_prob_smooth model.state_prob_smooth{n}];
end
[max_probs,post_states] = max(state_prob_smooth);
post_states = reshape(post_states',[no_obs,no_subjs]);
st = 1;
for subj = 1:no_subjs
    state_prob_smooth_subj = state_prob_smooth(:,st:st+no_obs-1);
    figure(subj)
    subplot(611)
    imagesc(state_prob_smooth_subj)
    title('State Probability Matrix')
    subplot(612)
    plot(task{subj},'k','linewidth',2)
    title('Actual Task')
    subplot(613)
    plot(est_states(:,subj),'r--','linewidth',2)
    title('Vitterbi Path')
    subplot(614)
    plot(post_states(:,subj),'g','linewidth',2)
    title('state(t)) = arg-max(state_prob_matrix(:,t)')
    counts_vitterbi = zeros(1,K);
    counts_post = zeros(1,K);
    for k = 1:K
        counts_vitterbi(k) = length(find(est_states(:,subj) == k));
        counts_post(k) = length(find(post_states(:,subj) == k));
    end
    subplot(615)
    stem(counts_vitterbi)
    title('Percent Occurence of each etate (Vitterbi Path)')
    subplot(616)
    stem(counts_post)
    title('Percent Occurence of each state (Posterior Prob) ')
    st = st + no_obs;
end
figure,
plot(model.F,'o-')
title('Log-Lower Bound of the best model')
figure
for k = 1:length(Model)
    subplot(length(Model),1,k)
    plot(Model{k}.F,'o-')
    title('Log-Lower Bound')
end
for k = 1:K
    ap = model.ap(k); bp = model.bp(:,:,k);
    est_cov(:,:,k) = bp/ap;
    invD = inv(diag(sqrt(diag(est_cov(:,:,k)))));
    pearson_corr(:,:,k) = invD*est_cov(:,:,k)*invD;
end
%estimeted state transition matrix A
Wa = model.Wa;
Wa = Wa';
H = size(Wa,1);
Aest = Wa./repmat(sum(Wa,1),H,1);% transition distribution p(h(t)|h(t-1))
%estimeted pi
Wpi = model.Wpi;
piest = Wpi./sum(Wpi);




