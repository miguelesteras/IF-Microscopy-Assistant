%% Phase 3. Create Orientation
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Determine the orientation of each cell based on the distance between
%   contour pixels and center of mass.
%   2. Compute center of mass for each individual cell, and distance from
%   center to contour pixels (radious)
%   3. Compute cell body and contour polar coordenates. 
%   ======================================================================
close all; clear; clc;

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                     
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    
    % index non-empty cells in cell sequences
    index = find(~cellfun(@isempty,cellSequences));
    
    % cell sequences index matches the new cell arrays created, containing
    % information computed for every non-empty cell.
    centerMass = cell(size(cellSequences));
    maxRadii = cell(size(cellSequences));
    rotationUp = cell(size(cellSequences));
    contourCoordinates = cell(size(cellSequences));
    cellCoordinates = cell(size(cellSequences));
    
    % for every non-empty cell in cell sequence; save pixel location 
    % (whole cell and contour) in polar system (to facilitate rotations) 
    % and normalized to center of mass at [0,0]. The rotation necessary to
    % make the contour pixel further from the center face North (theta=pi/2)
    % is calculate and will be used to normalize cell orientation before
    % training.
    for j = 1:size(index)
        selection = cellSequences{index(j)};
        image = false(metadata.imageSize);
        image(selection) = true;            
        stats = regionprops(image,'centroid','PixelList');  
        center = round(stats.Centroid);  
        centerMass{index(j)} = center;
        coordenates = stats.PixelList - ones(size(stats.PixelList))*diag(center);
        [theta, rho] = cart2pol(coordenates(:,1),coordenates(:,2)); 
        cellCoordinates{index(j)} = [wrapTo2Pi(theta) rho];
        [~, idx] = max(rho);
        rotation = pi/2 - wrapTo2Pi(theta(idx));
        maxRadii{index(j)} = [theta(idx) rho(idx)];
        rotationUp{index(j)} = rotation;
        SE = strel('disk',2);
        contour = bwmorph(imopen(image,SE),'remove');       
        stats2 = regionprops(contour,'PixelList');       
        carteCoordenates = stats2.PixelList - ones(size(stats2.PixelList))*diag(center);      
        [theta, rho] = cart2pol(carteCoordenates(:,1),carteCoordenates(:,2));
        contourCoordinates{index(j)} = [wrapTo2Pi(theta) rho];
    end
    
    maxRadious = ceil(max(cellfun(@(v) v(2), maxRadii(index))));    
    [metadata(:).maxRadious] = maxRadious;

    save(strcat(metadata.name,'_metadata.mat'),'metadata');
    save(strcat(metadata.name,'_centerMass.mat'),'centerMass');     
    save(strcat(metadata.name,'_maxRadii.mat'),'maxRadii');     
    save(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');     
    save(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     
    save(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates');     

    %clearvars -except files num_files i
end

%% Plots

% ps1 = polarscatter(wrapTo2Pi(theta),rho, 'filled');
% ps1.SizeData = 50;
% ps1.MarkerFaceAlpha = .5;
% hold on
% ps2 = polarscatter(wrapTo2Pi(theta+rotation),rho,'filled');
% ps2.SizeData = 50;
% ps2.MarkerFaceAlpha = .7;
% hold off
% lg = legend('Original','Normalized Rotation');
% lg.FontSize = 14;

