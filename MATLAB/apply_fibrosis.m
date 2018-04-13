function fibrotic = apply_fibrosis(pts, pattern)
%APPLY_FIBROSIS Summary of this function goes here
%   Detailed explanation goes here

fibrotic = apply_field(pts, pattern.oct_seed, pattern.region);
if isfield(pattern, 'dither')
    fibrotic = fibrotic & apply_field(pts, pattern.oct_seed, pattern.dither);
end
end

