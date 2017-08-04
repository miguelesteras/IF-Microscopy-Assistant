%% Phase 4.1. Transform cell silhouettes into Fourier descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell silhouette into Fourier descriptors.
%   Fourier descriptors as first described by Kuhl and Giardina in 
%   "Elliptic Fourier features of a closed contour" 
%   Computer Graphics and Image Processing 18:236-258 1982
%   This implementation has been created by David Thomas from the 
%   University of Melbourne, and can be found at:
%   https://uk.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors
%   Copyright (c) 2005, David Thomas
%   ======================================================================

%% Build data set of single cells

seqLength = 5;              % Number of frames in learning sequence (4 input + 1 target)
NoHarmonics = 25;           % Number of Harmonics (interger greater than 0)
NormSize = true;            % normalize size of object
NormOrientation = false;    % normalize orientation of object

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_contourCoordenates.mat'),'contourCoordenates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellSequences)),2);
    idx = find(noFrames >= seqLength);
    fourierInput = cell(100,4); 
    fourierTarget = cell(100,1);
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellSequences(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};
            m = 1;
            
            while m < seqLength                
                input = contourCoordenates{idx(j),idx2(k+m-1)};
                input(:,1) = input(:,1) + rotation;          % apply rotation
                [x,y] = pol2cart(input(:,1),input(:,2));
                contour = [x y];         
                descriptor = fEfourier(contour, NoHarmonics, NormSize, NormOrientation);
                fourierInput{count,m} = descriptor;
                m = m+1;
            end        

            target = contourCoordenates{idx(j),idx2(k+m-1)};
            target(:,1) = target(:,1) + rotation;          % apply rotation
            [x,y] = pol2cart(target(:,1),target(:,2));
            contour = [x y];         
            descriptor = fEfourier(contour, NoHarmonics, NormSize, NormOrientation);
            fourierTarget{count,1} = descriptor;
            count = count+1;
        end               
    end
    save(strcat(metadata.name,'_fourierInput.mat'),'fourierInput');     
    save(strcat(metadata.name,'_fourierTarget.mat'),'fourierTarget');     

    %clearvars -except files num_files k
end
