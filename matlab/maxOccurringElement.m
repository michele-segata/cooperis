function [v] = maxOccurringElement(array, random, seedn)

[counts, values] = groupcounts(array');
max_index = find(counts == max(counts));
if length(max_index) == 1
    v = values(max_index);
else
    if random
        rng(seedn);
        v = values(max_index(randi(length(max_index), 1)));
    else
        v = values(max_index(1));
    end
end