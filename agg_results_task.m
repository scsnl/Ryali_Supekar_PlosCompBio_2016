resultfmt = '/home/sryali/VB-HMM/Results/WM_HCP_1_20-%d.mat'
Model = {}

for ii=1:100
	foo = load(sprintf(resultfmt,ii));
	Model{ii} = foo.Model(:);
end

save('/home/sryali/VB-HMM/Results/WM_HCP_1_20-All.mat','Model')
