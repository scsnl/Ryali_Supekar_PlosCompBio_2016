function [fractional_occupancy, mean_life] = summary_stats_bystate(z,k)

T = length(z);
ix = find(z == k);
if ~isempty(ix)
    fractional_occupancy = length(ix)/T;
    counter = max(max(bwlabel(diag(z == k))));
    mean_life = length(ix)/counter;
else
    fractional_occupancy = 0;
    mean_life = 0;
end

