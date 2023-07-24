clear all
close all
clc

result_dir = '/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Scripts/multi_subject_combined/Results/';
std_noise = 0.25;
fname = strcat(result_dir,'sim_data_',num2str(std_noise),'.mat');
load(fname)
no_subjs = length(model.state_prob_smooth);
no_obs = size(model.state_prob_smooth{1},2);
K = size(model.Wa,1);
[est_states] =  est_states_vitterbi(data,model);
est_states = reshape(est_states',[no_obs,no_subjs]);
state_prob_smooth = [];
for n = 1:no_subjs
    state_prob_smooth = [state_prob_smooth model.state_prob_smooth{n}];
end
[max_probs,post_states] = max(state_prob_smooth);
post_states = reshape(post_states',[no_obs,no_subjs]);
st = 1;
figure,
plot(model.F,'o-')
for subj = 1:no_subjs
    state_prob_smooth_subj = state_prob_smooth(:,st:st+no_obs-1);
    figure(subj)
    subplot(611)
    imagesc(state_prob_smooth_subj)
    subplot(612)
    plot(task{subj},'k','linewidth',2)
    subplot(613)
    plot(est_states(:,subj),'r--','linewidth',2)
    subplot(614)
    plot(post_states(:,subj),'g','linewidth',2)
    counts_vitterbi = zeros(1,K);
    counts_post = zeros(1,K);
    for k = 1:K
        counts_vitterbi(k) = length(find(est_states(:,subj) == k));
        counts_post(k) = length(find(post_states(:,subj) == k));
    end
    subplot(615)
    stem(counts_vitterbi)
    subplot(616)
    stem(counts_post)
    st = st + no_obs;
end

for k = 1:K
    ap = model.ap(k); bp = model.bp(:,:,k);
    est_cov(:,:,k) = bp/ap;
end
%estimeted state transition matrix A
Wa = model.Wa;
Wa = Wa';
H = size(Wa,1);
Aest = Wa./repmat(sum(Wa,1),H,1);% transition distribution p(h(t)|h(t-1))
%estimeted pi
Wpi = model.Wpi;
piest = Wpi./sum(Wpi);
