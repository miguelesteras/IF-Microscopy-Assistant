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
    
    index = find(~cellfun(@isempty,cellSequences)); % non empty cells in cell sequences array
    centerMass = cell(size(cellSequences));
    maxRadii = cell(size(cellSequences));
    rotationUp = cell(size(cellSequences));

    for j = 1:size(index)
        selection = cellSequences{index(j)};    % select one cell
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'centroid');
        center = round(stats.Centroid);         % coordenates for center of mass
        SE = strel('disk',5);
        contour = bwmorph(imopen(image,SE),'remove');       
        stats = regionprops(contour,'PixelList');
        carteCoordenates = stats.PixelList - center;   % polar coordenates of cell contour (with origin in centroid)
        [theta, rho] = cart2pol(carteCoordenates(:,1),carteCoordenates(:,2));
        [radius, idx] = max(rho);
        rotation = -rad2deg(theta(idx)) + 90;
        centerMass{index(j)} = center;
        maxRadii{index(j)} = [theta(idx) radius];
        rotationUp{index(j)} = rotation;        
               
    end
    
    a = maxRadii{1,2};
    
    save(strcat(metadata.name,'_centerMass.mat'),'centerMass');     % save files on disk
    save(strcat(metadata.name,'_maxRadii.mat'),'maxRadii');     % save files on disk
    save(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');     % save files on disk

    clearvars -except files num_files k
end
