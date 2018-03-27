function f = find_level(x, threshold)
%FIND_LEVEL Percentage based scalar thresholding utility
%   F = FIND_LEVEL(X, THRESHOLD) Returns the value F such that
%   sum(X<=F)/length(X) is THRESHOLD
%
%   Uses a cumulative histogram with 1000 bins

    [bins, edges] = histcounts(x(:), 1000);
    nhist = cumsum(bins)/sum(bins);
    [~, k] = min(abs(nhist - threshold));
    f = edges(k+1);
end
    