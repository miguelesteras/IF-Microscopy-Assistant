%% Phase 3. Create Orientation
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
    contourCoordenates = cell(size(cellSequences));

    for j = 1:size(index)
        selection = cellSequences{index(j)};    % select one cell
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'centroid');
        center = round(stats.Centroid);         % coordenates for center of mass
        SE = strel('disk',2);
        contour = bwmorph(imopen(image,SE),'remove');       
        stats = regionprops(contour,'PixelList');
        carteCoordenates = stats.PixelList - center;        % cartesian coordenates of cell contour (with origin in centroid)
        contourCoordenates{index(j)} = carteCoordenates;
        [theta, rho] = cart2pol(carteCoordenates(:,1),carteCoordenates(:,2));   % polar coordenates
        [~, idx] = max(rho);
        rotation = pi/2 - theta(idx);
        centerMass{index(j)} = center;
        maxRadii{index(j)} = [theta(idx) rho(idx)];
        rotationUp{index(j)} = rotation;
        contourCoordenates{index(j)} = [theta rho]; % polar coordenates

    end
    
    maxRadious = ceil(max(cellfun(@(v) v(2), maxRadii(index))));    
    [metadata(:).maxRadioius] = maxRadious;

    save(strcat(metadata.name,'_metadata.mat'),'metadata');
    save(strcat(metadata.name,'_centerMass.mat'),'centerMass');     
    save(strcat(metadata.name,'_maxRadii.mat'),'maxRadii');     
    save(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');     
    save(strcat(metadata.name,'_contourCoordenates.mat'),'contourCoordenates');     

    clearvars -except files num_files k
end
