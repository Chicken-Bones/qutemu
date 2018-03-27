function noise = apply_atria(pts, table, unit, harmonics, roughness, fill)
%APPLY_ATRIA Applies an atrial perlin noise field to a set of a points
%   noise = APPLY_ATRIA(PTS, TABLE, UNIT, HARMONICS, ROUGHNESS, FILL)
%   PTS a N-by-3 matrix of points to sample the noise distribution at
%   TABLE a struct array used as seed returned from gen_octave()
%   UNIT the wavelength of the base spatial frequency, the grid size of
%   first harmonic
%   HARMONICS the number of noise harmonics to use (1 to length(table))
%   ROUGHNESS the relative amplitude of each successive harmonic (1.0 is
%   all harmonics have equal weight, 0.0 is just the first harmonic)
%   FILL a ratio for thresholding the noise to a logical array
%
%   noise = APPLY_ATRIA(PTS, TABLE, CFG) packs the arguments above into a
%   structure CFG

    if isa(unit, 'struct')
        cfg = unit;
        unit = cfg.unit;
        harmonics = cfg.harmonics;
        roughness = cfg.roughness;
        fill = cfg.fill;
    end

    noise = octave(pts / unit, table(1:harmonics), roughness);
    if (fill > 0)
        noise = noise <= find_level(noise, fill);
    end
end