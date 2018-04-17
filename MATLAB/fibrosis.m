folder = fileparts(mfilename('fullpath'));
folder = fullfile(folder, '../data/fibrosis');

nodes = LoadNodeCoords;
ele = LoadEleIndices;
elepos = permute(reshape(nodes(ele, :)', 3, [], size(ele, 2)), [2 1 3]);

load(fullfile(folder, 'pattern')); %fibrotic_pattern
fills = [0.6 0.7 0.8, 0.9, 0.95, 1];
for fill = fills
    fprintf('Applying %d%%\n', fill*100);
    pattern = fibrotic_pattern;
    pattern.dither.fill = fill;
    fib = apply_fibrosis(elepos, pattern);
    
    h5path = fullfile(folder, ['ptn0fill', char(string(fill*100)), '.h5']);
    h5create(h5path,'/Conductivity', size(fib))
    h5write(h5path, '/Conductivity', single(1 - fib))
end