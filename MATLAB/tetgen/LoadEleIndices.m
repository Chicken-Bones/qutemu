function [ele, meta] = LoadEleIndices(path)
%LoadEleIndices Tetgen .ele file loader
%   LoadEleIndices(path) Loads the vertex data in a tetgen ele file into a
%   MATLAB array. Vertices are returned as MATLAB indices (starting at 1)
%   into a corresponding points table
%
%   LoadEleIndices() Assumes the path is ../../data/Karlsruhe/mesh.ele

if nargin == 0
    folder = fileparts(mfilename('fullpath'));
    path = fullfile(folder, '../../data/Karlsruhe/mesh.ele');
end

file_obj=fopen(path);
read_data=cell2mat(textscan(file_obj,'','headerlines',1,'delimiter',' ','collectoutput',1));
fclose(file_obj);

ele = read_data(:,2:5)+1;
meta = read_data(:, 6:end);

end

