function noise = octave(pts, tables, persistence, twoD)
%OCTAVE 3D Octave Perlin Noise Generator
%   OCTAVE(PTS, TABLES, PERSISTANCE) Calculates octave noise at the
%   coordinates specified in PTS using permutation TABLES provided by
%   seed_octave(n)
%   
%   Octave noise is a scaled sum of increasing resolution perlin noise The
%   number of octaves is length(TABLES) Each octave table contains an
%   offset point in the unit cube which is added to the points to remove
%   grid alignment Points for octave I are multiplied by 2^(i-1), halving
%   the grid size each octave Noise for octave I is scaled by
%   PERSISTANCE^(i-1) so higher resolutions (details) have less weight
%
%   OCTAVE(PTS, TABLES) Uses the default value 0.5 for PERSISTANCE
%
%   OCTAVE(PTS, TABLES, PERSISTANCE, true) Sets the z-offset in the
%   permutation tables to 0, so that 2D noise can be generated with integer
%   z coordinates

    if nargin < 3
        persistence = 0.5;
    end
    
    n = length(tables);
    if nargin >= 4 && twoD %need integer z coordinate for 2d noise to have reasonable spatial properties
        for i=1:n
            tables(i).offset(3) = 0;
        end
    end
    
    noise = zeros(size(pts, 1), 1);
    for i = 1:n
        tbl = tables(i);
        p = i-1;
        noise = noise + perlin(pts * 2^p + tbl.offset, tbl.table) * persistence^p;
    end
    
    max = sum(persistence.^(0:n-1));
    noise = noise / max;
end