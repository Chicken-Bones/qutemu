folder = fileparts(mfilename('fullpath'));
folder = fullfile(folder, '../data/fibrosis');
if exist(folder, 'dir') ~= 7
    mkdir(folder)
end


nodes = LoadNodeCoords;
ele = LoadEleIndices;
elepos = mean(permute(reshape(nodes(ele, :)', 3, [], size(ele, 2)), [2 1 3]), 3);

sa_center = [8.02126 7.81824 40.1448];
sa_radius = 1.116;
sa_region = sqrt(sum((elepos - sa_center).^2, 2)) < sa_radius;

if exist('fibrotic_pattern', 'var') ~= 1
    load(fullfile(folder, 'pattern')); %fibrotic_pattern
end
fills = [0.6 0.7 0.8, 0.9, 0.95, 1];
for fill = fills
    fprintf('Applying %d%%\n', fill*100);
    pattern = fibrotic_pattern;
    pattern.dither.fill = fill;
    fib = apply_fibrosis(elepos, pattern);
    fib(sa_region) = 0;
    
    h5path = fullfile(folder, ['ptn0fill', char(string(fill*100)), '.h5']);
    if exist(h5path, 'file') == 2
        delete(h5path)
    end
    h5create(h5path,'/Conductivity', size(fib))
    h5write(h5path, '/Conductivity', single(1 - fib))
end
save(fullfile(folder, 'pattern'), 'fibrotic_pattern')