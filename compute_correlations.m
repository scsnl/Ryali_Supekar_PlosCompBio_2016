function [correlations,zcorrelations] = compute_correlations(data_subj,Labels_subj,label)

%Compute correlations for each subject at each state
no_subjects = length(data_subj);
zcorrelations = [];
correlations = [];
cnt = 1;
for subj = 1:no_subjects
    %X = data_subj(:,:,subj);
    X = data_subj{subj};
    %ix = find(Labels_subj(:,subj) == label);
    ix = find(Labels_subj{subj} == label);
    if ~isempty(ix)
        Xlabel = X(:,ix);
        fc = corr(Xlabel');
        fcz = log((1+fc)./(1-fc));
        zcorrelations(:,:,cnt) = fcz;
        correlations(:,:,cnt) = fc;
        cnt = cnt + 1;
    end
end
