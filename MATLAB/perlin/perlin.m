function [noise] = perlin(pts, P)
%PERLIN 3D Perlin Noise Generator
%   PERLIN(PTS, TABLE) Calculates perlin noise at the coordinates specified
%   in PTS using a permutation TABLE provided by seed_perlin()
%   
%   The implementation is vectorised, but still slow on large datasets due
%   to large memory footprints and slow bitwise operations in MATLAB
%
%   This function generates 3D noise, if 2D noise is desired, use an
%   integer coordinate for Z
%
%   Based on https://mrl.nyu.edu/~perlin/paper445.pdf 
%       (K Perlin. 2002. Improving noise. In SIGGRAPH, 681–682.)
%   Hash from https://www.cs.utah.edu/~aek/research/noise.pdf
%       (A Kensler, A Knoll, P Shirley. 2008. Better Gradient Noise. SCI Institute Technical Report No. UUSCI-2008-001)

% Use 256 evenly distributed vectors on the unit sphere for gradients
% Original implementations just use 16
sphere = sphere_fibonacci_grid_points(256);
swizzle = ((dec2bin(0:7)-'0')*3 + [1 2 3]);

% integer and fractional parts of grid
boxmin = floor(pts);
boxmax = ceil(pts);
t = pts - boxmin;

% select gradient vectors at box corners
AB = [perlin_permute(boxmin) perlin_permute(boxmax)];
K = reshape(AB(:, swizzle), [], 8, 3);
k = bitxor(bitxor(K(:, :, 1), K(:, :, 2)), K(:, :, 3))+1;
grads = permute(reshape(sphere(k', :), 8, [], 3), [2, 1, 3]);

% vectors from pt to each corner
dcnrs = [t t-1];
dcnrs = reshape(dcnrs(:, swizzle), [], 8, 3);

% noise factors, gradient vector dot products
noise = sum(grads .* dcnrs, 3);

% apply fade
t = t .* t .* t .* (t .* (t * 6 - 15) + 10);

% lerp between the 8 sample points
noise = noise(:, 1:4).*(1-t(:, 1)) + noise(:, 5:8).*t(:, 1);
noise = noise(:, 1:2).*(1-t(:, 2)) + noise(:, 3:4).*t(:, 2);
noise = noise(:, 1  ).*(1-t(:, 3)) + noise(:, 2  ).*t(:, 3);

% Lookup each coordinate of i in the permutation tables
function [p] = perlin_permute(i)
    i = mod(i, 256) + 1;
    p = [P(i(:, 1), 1) P(i(:, 2), 2) P(i(:, 3), 3)];
end

end