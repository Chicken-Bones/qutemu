function [params] = seed_octave(n)
%SEED_OCTAVE Octave Noise Seed Generator
%   SEED_OCTAVE(N) Generates an vector of N structs containing a perlin
%   table and an offset vector in the unit cube

    tables = cell(n, 1);
    offsets = cell(n, 1);
    for i = 1:n
        tables{i} = seed_perlin();
        offsets{i} = rand(1, 3);
    end
    params = struct('table', tables, 'offset', offsets);
end