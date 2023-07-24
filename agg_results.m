resultfmt = '/home/sryali/VB-HMM/Results/stanford_dataset_grp2-%d.mat'
Model = {}

for ii=1:100
	foo = load(sprintf(resultfmt,ii));
	Model{ii} = foo.Model(:);
end

save('/home/sryali/VB-HMM/Results/stanford_dataset_grp2-All.mat','Model')
%save('/home/ksupekar/stanford_dataset_grp2-All.mat','Model')
