function main_vbhmm_HCP_replication_session2(ss)

warning('off')
addpath('/home/sryali/VB-HMM/Scripts')

K = 25; % # of states
maxIter = 100; % maximum number of allowed iterations
tol = 10^-3;   % tolerance for convergence
init_option = 3; %1-random, 2-task,3-kmeans

% Please set an output directory and unique filename
saveDir = '/home/sryali/VB-HMM/Results/';
result_fname = sprintf('Adults_6ROIs_25Clusters_dataset1_session2-%d.mat',ss);

% Input data filename
inFile = '/home/sryali/VB-HMM/Data/data_set1_session2.mat';

% Get the  Data
task = [];
load(inFile)

% Setup parcluster -- create unique job directory
cwd = pwd;
mkdir(fullfile(cwd,'temp'));
jobDir = fullfile(cwd,'temp',sprintf('job%d',ss));
if exist(jobDir,'dir')
	rmdir(jobDir,'s');
end
mkdir(jobDir);
pc = parcluster()
pc.NumWorkers = str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
pc.JobStorageLocation = jobDir;
delete(pc.Jobs)

%% We delay for 0-5s because simultaneous starts of the parcluster cause problems
% 	even when the job directory is unique
rng(ss);
delay = 5*rand(1); %0-5s delay
pause(delay)
rng default; %reset the rng

matlabpool(pc)
[Model]= vb_cont_hmm_multi_subj_1rep(data,K,maxIter,tol,task,init_option,ss);
save([saveDir result_fname], 'Model')
delete(gcp('nocreate'))
rmdir(jobDir)

