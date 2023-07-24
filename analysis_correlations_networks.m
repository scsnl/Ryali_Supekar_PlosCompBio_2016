clear all
close all
clc

addpath(genpath('/mnt/musk2/home/fmri/fmrihome/scripts/GCCA_toolbox_sep21'))
addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts'));
addpath(genpath('/mnt/mapricot/musk2/home/sryali/Work/toolboxes/VB_HMM/vbhmm'))
addpath('/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Scripts/common_scripts')
addpath(genpath('/home/tianwenc/Toolbox/BCT/BCT_04_05_2014'));
data_dir = '/mnt/mandarin2/Public_Data/HCP/Stats/TS_for_TRSBN/Smith/';
prefix = 'NonInterestRegIC_eigen1_ts_rfMRI_REST1_LR_';
saveDir = '/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Results/HCP/';
result_fname = 'Adults_6ROIs_25Clusters_lambda_10p3_rep2.mat';


load([saveDir result_fname])
subjects = ReadList('/mnt/mandarin1/Public_Data/HCP/Data/subjectslist.txt');
roi_names =  {'ACC','lAI','rAI','lDLPFC','lPPC','rDLPFC','rPPC','preCue','VMPFC'};
rois = [1,3,6,7,8,9];
roi_names = roi_names(rois);
data = [];


data = [];
pval = 10^-3;
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
K = size(model.Wa,1);
counts_post = zeros(1,K);
for k = 1:K
    counts_post(k) = length(find(post_states == k));
end
counts_post = counts_post/sum(counts_post);
[percent_dominant dominant_states] = sort(counts_post,'descend');
figure
bar(counts_post*100)
%[M,T,S] = size(data_subj);
%Labels_subj = reshape(post_states,T,S);
st = 1;
for s = 1:length(data_states)
    T = size(data_states{s},2);
    Labels_subj{s} = post_states(st:st+T-1);
    st = st + T;
end

for k = 1:3
    label = dominant_states(k);
    [correlations_subj,zcorrelations_subj] = compute_correlations(data_states,Labels_subj,label);
    [H,P,CI,STATS] = ttest(zcorrelations_subj,0,pval,'right',3);
    figure(2)
    subplot(1,3,k)
    cca_plotcausality(H,roi_names,5);
     tstat = STATS.tstat;
     ix = isinf(tstat); tstat(ix) = 1;
     [est_network,clust_mtx] = clusters_community_detection(tstat);
%    avg_cov = mean(correlations_subj,3);
     figure(3)
    subplot(1,3,k)
    cca_plotcausality(est_network,roi_names,5);
end
