function [P] = seed_perlin()
%SEED_PERLIN Perlin Noise Seed Generator
%   SEED_PERLIN() Generates a permutation table 0-255 for x,y,z
    P = [randperm(256)' randperm(256)' randperm(256)']-1;
end