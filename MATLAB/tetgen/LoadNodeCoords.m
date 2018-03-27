function pts = LoadNodeCoords(path)
%LoadNodeCoords Tetgen .ele file loader
%   LoadNodeCoords(path) Loads the vertex data in a tetgen node file into a
%   MATLAB array [X1 Y1 Z1; X2 Y2 Z2; ...
%
%   LoadNodeCoords() Assumes the path is ../../data/Karlsruhe/mesh.node

if nargin == 0
    folder = fileparts(mfilename('fullpath'));
    path = fullfile(folder, '../../data/Karlsruhe/mesh.node');
end

file_obj=fopen(path);
read_data=cell2mat(textscan(file_obj,'','headerlines',1,'delimiter',' ','collectoutput',1));
fclose(file_obj);

pts = read_data(:,2:4);

end

