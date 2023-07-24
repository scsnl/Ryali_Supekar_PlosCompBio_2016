clear all
close all
clc

addpath(genpath('/mnt/musk2/home/fmri/fmrihome/scripts/GCCA_toolbox_sep21'))
addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts'));
addpath(genpath('/mnt/mapricot/musk2/home/sryali/Work/toolboxes/VB_HMM/vbhmm'))
addpath('/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Scripts/common_scripts')
addpath(genpath('/home/tianwenc/Toolbox/BCT/BCT_04_05_2014'));


data_dir = '/mnt/mandarin1/Public_Data/HCP/Data_set2_20ss/TS_for_TRSBN/ICA_REST1_NonInterestRegIC/';
prefix = 'regWMCSF_ms_eigen1_ts_rfMRI_REST1_LR_';
saveDir = '/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Results/HCP/';
result_fname = 'Adults_9ROIs_25Clusters_lambda_10p3_replication_data.mat';

load([saveDir result_fname])
subjects = ReadList('/mnt/mandarin1/Public_Data/HCP/Data_set2_20ss/subjectslist.txt');
roi_names =  {'ACC','lAI','rAI','lDLPFC','lPPC','rDLPFC','rPPC','preCue','VMPFC'};
%rois = [1,3,6,7,8,9];
rois = 1:9;
roi_names = roi_names(rois);
data = [];
%%%%%%%%%%% Get the  Data %%%%%%%%%%%%%%%
for subj = 1:length(subjects)
    load([data_dir prefix subjects{subj} '.mat']);
    X = roi_data.timeseries;
    X = X(rois,:);
    for k = 1:size(X,1)
        x = X(k,:);
        X(k,:) = (x-mean(x))/std(x);
    end
    data_subj(:,:,subj) = X; %for Vitterbi states
    data_states{subj} = X;
end
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
    %Partial Correlation
    inv_est_cov(:,:,k) = inv(est_cov(:,:,k));
    invD = inv(diag(sqrt(diag(inv_est_cov(:,:,k)))));
    partial_corr(:,:,k) = -invD*inv_est_cov(:,:,k)*invD;
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
for k = 1:6
    [est_network1,clust_mtx1] = clusters_community_detection((pearson_corr(:,:,dominant_states(k))));
   [est_network2,clust_mtx2] = clusters_community_detection((partial_corr(:,:,dominant_states(k))));
   figure(5)
    subplot(2,3,k)
    cca_plotcausality(est_network1,roi_names,5);
    figure(6)
    subplot(2,3,k)
    cca_plotcausality(est_network2,roi_names,5);
end
static_corr = corr(cell2mat(data_states)');
invD = inv(diag(sqrt(diag(static_corr))));
inv_cov = inv(cov(cell2mat(data_states)'));
invD = inv(diag(sqrt(diag(inv_cov))));
partial_corr = invD*inv_cov*invD;
[est_network1,clust_mtx1] = clusters_community_detection((static_corr));
[est_network2,clust_mtx2] = clusters_community_detection((partial_corr));
figure(7)
subplot(211)
cca_plotcausality(est_network1,roi_names,5);
subplot(212)
cca_plotcausality(est_network2,roi_names,5);

