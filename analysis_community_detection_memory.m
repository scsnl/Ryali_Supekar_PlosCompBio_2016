clear all
close all
clc

addpath(genpath('/mnt/musk2/home/fmri/fmrihome/scripts/GCCA_toolbox_sep21'))
addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts'));
addpath(genpath('/mnt/mapricot/musk2/home/sryali/Work/toolboxes/VB_HMM/vbhmm'))
addpath('/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Scripts/common_scripts')
addpath(genpath('/mnt/mapricot/musk2/home/sryali/Work/toolboxes/VB_HMM/vbhmm'))
addpath(genpath('/home/tianwenc/Toolbox/BCT/BCT_04_05_2014'));

%data_dir = '/mnt/mabloo1/apricot1_share6/memory_consolidation/dynamic_states/ROIs/TimeSeries/';
data_dir = '/mnt/apricot1_share6/memory_consolidation/dynamic_states/Results/TimeSeries_RS2/';
saveDir = '/mnt/apricot1_share6/memory_consolidation/dynamic_states/Results/';
result_fname = 'children_post_encoding_dynamic_states_memory_RS2.mat';
load([saveDir result_fname])
subjects = ReadList('/mnt/mandarin1/Public_Data/HCP/Data/subjectslist.txt');
data = [];
%%%%%%%%%%% Get the  Data %%%%%%%%%%%%%%%
subjects = ReadList([data_dir,'subjectlist.txt']);
%%%%%%%%%%% Get the  Data %%%%%%%%%%%%%%%
for subj = 1:length(subjects)
       load(subjects{subj});
        X = all_roi_ts';
        for k = 1:size(X,1)
            x = X(k,:);
            X(k,:) = (x-mean(x))/std(x);
        end
        data_subj(:,:,subj) = X(:,1:170); %for Vitterbi states
        data_states{subj} = X;
end
roi_names = {'lHIPP', 'lIFG', 'lPPC', 'lPPA', 'rHIPP', 'rIFG', 'rPPA','rPPC'};
figure(1)
for k = 1:length(Model)
    %subplot(length(Model),1,k)
    plot(Model{k}.F + 5*randn,'o-')
    hold on
    title('Log-Lower Bound')
end
for repetition = 1:length(Model)
    F = Model{repetition}.F;
    if repetition == 1
        model = Model{repetition};
    elseif (max(model.F) < max(F))
        model = Model{repetition};
    end
end
state_prob_smooth = [];
Task = [];
for n = 1:length(data_states)
    state_prob_smooth = [state_prob_smooth model.state_prob_smooth{n}];
    %    % Task = [Task task{n}];
end
[est_states] =  est_states_vitterbi(data_states,model);
post_states = est_states;

%[max_probs,post_states] = max(state_prob_smooth);

figure(2)
subplot(311)
imagesc(state_prob_smooth)
subplot(312)
plot(post_states,'linewidth',2)
subplot(313)
plot(model.F,'o-','linewidth',2)

K = size(model.Wa,1);
counts_post = zeros(1,K);
for k = 1:K
    counts_post(k) = length(find(post_states == k));
end
counts_post = counts_post/sum(counts_post);
[percent_dominant dominant_states] = sort(counts_post,'descend');
figure(3)
bar(counts_post*100)
%Estimates of Covariance
for k = 1:K
    ap = model.ap(k); bp = model.bp(:,:,k);
    est_cov(:,:,k) = bp/ap;
    invD = inv(diag(sqrt(diag(est_cov(:,:,k)))));
    pearson_corr(:,:,k) = invD*est_cov(:,:,k)*invD;
    inv_est_cov(:,:,k) = inv(est_cov(:,:,k));
    invD = inv(diag(sqrt(diag(inv_est_cov(:,:,k)))));
    partial_corr(:,:,k) = invD*inv_est_cov(:,:,k)*invD;
end
%estimeted state transition matrix A
Wa = model.Wa;
Wa = Wa';
H = size(Wa,1);
Aest = Wa./repmat(sum(Wa,1),H,1);% transition distribution p(h(t)|h(t-1))
%estimeted pi
Wpi = model.Wpi;
piest = Wpi./sum(Wpi);
figure(4)
for k = 1:15
    subplot(3,5,k)
    cca_plotcausality(abs(pearson_corr(:,:,dominant_states(k))) > 0.2,roi_names,5);
end
figure(5)
for k = 1:15
    %[est_network,clust_mtx] = clusters_community_detection((pearson_corr(:,:,dominant_states(k))));
    [est_network,clust_mtx] = clusters_community_detection((partial_corr(:,:,dominant_states(k))));
    subplot(3,5,k)
    cca_plotcausality(est_network,roi_names,5);
end
static_corr = inv(corr(cell2mat(data_states)'));
inv_static_corr = inv(static_corr);
invD = inv(diag(sqrt(diag(inv_static_corr))));
static_partial_corr = invD*inv_static_corr*invD;
%[est_network,clust_mtx] = clusters_community_detection((static_corr));
[est_network,clust_mtx] = clusters_community_detection((static_partial_corr));
figure(6)
cca_plotcausality(est_network,roi_names,5);
