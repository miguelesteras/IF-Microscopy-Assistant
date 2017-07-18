%% Phase 3. Create Orientation (North point)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Determine the orientation of each cell based on the distance between
%   contour pixels and center of mass. 
%   ======================================================================

%% Build data set of single cells

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                     % load files from disk
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');

    northPoints = cell(size(cellSequences));        % initialize array for north point coordenates for full view singl
    index = find(~cellfun(@isempty,cellSequences)); % non empty cells in cell sequences array
    
    for j = 1:size(index)
        selection = cellSequences{index(j)};    % select one cell
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'centroid');
        center = round(stats.Centroid);         % coordenates for center of mass
        SE = strel('disk',5);
        contour = bwmorph(imopen(image,SE),'remove');       
        stats = regionprops(contour,'PixelList');
        coordenates = stats.PixelList;          % coordenates of cell contour 
        L2norm = sqrt((coordenates(:,1) - center(1)).^2 + (coordenates(:,2) - center(2)).^2);   % Euclidian distance from contour pixels to center of mass
        [~,idx] = max(L2norm);
        point = coordenates(idx,:);             % coordenates of contour pixel further from center of mass
        northPoints{index(j)} = point;           
    end
    
    save(strcat(metadata.name,'_northPoints.mat'),'northPoints');     % save files on disk
    
end
