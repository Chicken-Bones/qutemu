folder = fileparts(mfilename('fullpath'));
folder = fullfile(folder, '../data/Karlsruhe/');

[pts, pmeta] = LoadNodeCoords(fullfile(folder, 'mesh_original.node'));
n = size(pts, 1);
ele = LoadEleIndices();
connmap = cell(n, 1);
for i = 1:size(ele, 1)
    e = ele(i, :);
    connmap{e(1)} = [connmap{e(1)} e([2 3 4])];
    connmap{e(2)} = [connmap{e(2)} e([1 3 4])];
    connmap{e(3)} = [connmap{e(3)} e([1 2 4])];
    connmap{e(4)} = [connmap{e(4)} e([1 2 3])];
end
connmap = cellfun(@unique, connmap, 'UniformOutput', false);
nconn = cellfun(@length, connmap);

stim = pmeta(:, 2);
for site = 1:max(stim)
    isstim = stim == site;
    stimconn = cellfun(@(a) sum(isstim(a)), connmap);
    fprintf('start:%d\n', sum(isstim))

    % nodes with more than 70% connection to other stim nodes are in
    while 1
        nstim = sum(isstim);
        isstim = isstim | stimconn > nconn * 0.7;
        if sum(isstim) == nstim
            break
        end
        stimconn = cellfun(@(a) sum(isstim(a)), connmap);
        fprintf('dilate:%d\n', sum(isstim))
    end

    % nodes with less than 20% connection to other stim nodes are out
    while 1
        nstim = sum(isstim);
        isstim = isstim & stimconn > nconn * 0.2;
        if sum(isstim) == nstim
            break
        end
        stimconn = cellfun(@(a) sum(isstim(a)), connmap);
        fprintf('erode:%d\n', sum(isstim))
    end
    
    stim(stim == site) = 0;
    stim(isstim) = site;
end

h5path = fullfile(folder, 'stim.h5');
delete(h5path)
h5create(h5path,'/Stim', size(stim))
h5write(h5path, '/Stim', int8(stim))

pmeta(:, 2) = stim;
file = fopen(fullfile(folder, 'mesh.node'),'w');
fprintf(file,'%d\t%d\t%d\t0\n', n, size(pts, 2), size(pmeta, 2));
fmat = [(0:n-1)', pts, pmeta];
fprintf(file,'%d %f %f %f %d %d\n',fmat'); %hardcode format string for now
fclose(file);