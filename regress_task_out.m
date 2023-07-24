function yres = regress_task_out(y,stim_design,TR)

%Regress the out the task

[M,N] = size(y);
yres = zeros(N,M);
X(:,1) = ones(N,1);
hrf = spm_hrf(TR);
regressor = conv(hrf,stim_design)';
X(:,2) = regressor(1:N);
y = y';
for m = 1:M
    yres(:,m) = y(:,m) - X*(X\y(:,m));
    yres(:,m) = yres(:,m) - mean(yres(:,m))/std(yres(:,m));
end

