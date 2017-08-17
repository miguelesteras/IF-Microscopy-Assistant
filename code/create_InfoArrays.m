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


files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                     % load files from disk
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    
    % index non-empty cells in cell sequences
    index = find(~cellfun(@isempty,cellSequences));
    
    % cell sequences index matches the new cell arrays created, containing
    % information computed for every non-empty cell.
    centerMass = cell(size(cellSequences));
    maxRadii = cell(size(cellSequences));
    rotationUp = cell(size(cellSequences));
    contourCoordenates = cell(size(cellSequences));
    cellCoordenates = cell(size(cellSequences));
    
    % for every non-empty cell in cell sequence; save pixel location 
    % (whole cell and contour) in polar system (to facilitate rotations) 
    % and normalized to center of mass at [0,0]. The rotation necessary to
    % make the contour pixel further from the center face up (theta=pi/2)
    % is calculate and will be used to normalize cell orientation before
    % training.
    for j = 1:size(index)
        selection = cellSequences{index(j)};
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'centroid','PixelList');  
        center = round(stats.Centroid);            
        coordenates = stats.PixelList - center;
        [theta, rho] = cart2pol(coordenates(:,1),coordenates(:,2)); 
        cellCoordenates{index(j)} = [theta rho];        
        SE = strel('disk',2);
        contour = bwmorph(imopen(image,SE),'remove');       
        stats2 = regionprops(contour,'PixelList');       
        carteCoordenates = stats2.PixelList - center;      
        contourCoordenates{index(j)} = carteCoordenates;
        [theta, rho] = cart2pol(carteCoordenates(:,1),carteCoordenates(:,2));
        [~, idx] = max(rho);
        rotation = pi/2 - theta(idx);
        centerMass{index(j)} = center;
        maxRadii{index(j)} = [theta(idx) rho(idx)];
        rotationUp{index(j)} = rotation;
        contourCoordenates{index(j)} = [theta rho];

    end
    
    maxRadious = ceil(max(cellfun(@(v) v(2), maxRadii(index))));    
    [metadata(:).maxRadious] = maxRadious;

    save(strcat(metadata.name,'_metadata.mat'),'metadata');
    save(strcat(metadata.name,'_centerMass.mat'),'centerMass');     
    save(strcat(metadata.name,'_maxRadii.mat'),'maxRadii');     
    save(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');     
    save(strcat(metadata.name,'_contourCoordenates.mat'),'contourCoordenates');     
    save(strcat(metadata.name,'_cellCoordenates.mat'),'cellCoordenates');     

    clearvars -except files num_files i
end
