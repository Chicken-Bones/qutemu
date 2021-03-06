function [data, t, nodemap] = load_results(path)
%LOAD_RESULTS Summary of this function goes here
%   Detailed explanation goes here

h5path = fullfile(path, 'results.h5');
data = squeeze(h5read(h5path, '/Data'));
t = h5read(h5path, '/Data_Unlimited');
nodemap = (1:size(data, 1))-1;
if ~h5readatt(h5path, '/Data', 'IsDataComplete')
    nodemap = h5readatt(h5path, '/Data', 'NodeMap');
end

ppath = fullfile(path, 'permutation.txt');
if exist(ppath, 'file')
    file_obj=fopen(ppath);
    perm=cell2mat(textscan(file_obj,'','headerlines',1,'delimiter',' ','collectoutput',1));
    fclose(file_obj);
    
    if ~isempty(nodemap)
        iperm = zeros(size(perm, 1), 1);
        iperm(perm(:, 2)+1) = 1:length(iperm);
        nodemap = iperm(nodemap+1)-1;
        [nodemap, k] = sort(nodemap);
        data = data(k, :);
    else
        data = data(perm(:, 2)+1, :);
    end
end

