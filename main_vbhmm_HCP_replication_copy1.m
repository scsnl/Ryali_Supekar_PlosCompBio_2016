clear all
close all
clc
warning('off')
addpath(genpath('/mnt/mapricot/musk2/home/sryali/Work/toolboxes/VB_HMM/vbhmm'))
addpath('/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Scripts/common_scripts')
addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts'));

K = 25; % # of states
maxIter = 100; % maximum number of allowed iterations
tol = 10^-3;   % tolerance for convergence
no_repetitions = 100; % Run vbhmm 100 times and choose the best model
init_option = 3; %1-random, 2-task,3-kmeans

data_dir = '/mnt/mandarin1/Public_Data/HCP/Data_set3_9ss/TS_for_TRSBN/ICA_REST1_NonInterestRegIC/';
prefix = 'ms_eigen1_ts_rfMRI_REST1_LR_';
saveDir = '/mnt/mapricot/musk2/home/sryali/Work/switchingFC/VB-HMM/Results/HCP/';
result_fname = 'Adults_9ROIs_25Clusters_lambda_10p3_replication_data_no_reg_WMCSF_9commonsubjs.mat';
%%%%%%%%%%% Get the  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjects = ReadList('/mnt/mandarin1/Public_Data/HCP/Data_set3_9ss/subjectslist_common.txt');
roi_names =  {'ACC','lAI','rAI','lDLPFC','lPPC','rDLPFC','rPPC','preCue','VMPFC'};
%rois = [1,3,6,7,8,9];
rois = 1:9;
%%%%%%%%%%% Get the  Data %%%%%%%%%%%%%%%
for subj = 1:length(subjects)
    load([data_dir prefix subjects{subj} '.mat']);
    X = roi_data.timeseries;
    X = X(rois,:);
    for k = 1:size(X,1)
        x = X(k,:);
        X(k,:) = (x-mean(x))/std(x);
    end
    data{subj} = X;
    %subplot(4,5,subj)
    %plot(X')
end
task = [];
% VB-HMM
matlabpool local 8
[Model]= vb_cont_hmm_multi_subj_test(data,K,maxIter,tol,no_repetitions,task,init_option);
matlabpool close
save([saveDir result_fname], 'Model')



