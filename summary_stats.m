function [fractional_occupancy, mean_life,Counters] = summary_stats(z,K)

T = length(z);
fractional_occupancy = zeros(1,K);
mean_life = zeros(1,K);
Counters = zeros(1,K);
for k = 1:K
    ix = find(z == k);
    if ~isempty(ix)
        fractional_occupancy(k) = length(ix)/T;
    end
%     flag = 1;
%     counter = 0;
%     curri = 1;
%     while curri < length(ix)
%         if (ix(curri+1)  - ix(curri) ~= 1)
%             counter = counter + 1;
%         end
%            curri = curri+1;
%     end
    counter = max(max(bwlabel(diag(z == k))));
    Counters(k) = counter;
    mean_life(k) = length(ix)/counter;
end
