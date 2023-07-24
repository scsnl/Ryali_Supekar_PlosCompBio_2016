function [x,task] = sim_corr_data(stim_design,std_noise)

% Simulate 2-node network for a Block Design
C1 = [1 0;0 1];
mu1 = [0 0];
C2 = [1 -1;-1 1];
mu2 = [0 0];
x = zeros(2,length(stim_design));
for t = 1:length(stim_design)
    if stim_design(t) == 0
        x(:,t) = mvnrnd(mu1,C1);
    else
        x(:,t) = mvnrnd(mu2,C2);
    end
end
x = x +randn(size(x)).*std_noise;
task = zeros(1,length(stim_design));
task(stim_design == 0) = 1;
task(stim_design == 1) = 2;

%smooth the data with a moving average filter
% a = 1; b = ones(1,9).*1/9;
% for m = 1:size(x,1)
%     x(m,:) = filter(b,a,x(m,:));
% end

