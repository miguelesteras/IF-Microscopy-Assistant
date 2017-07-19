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
%   ======================================================================

%% Build data set of single cells

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                     % load files from disk
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    load(strcat(metadata.name,'_northPoints.mat'),'northPoints');

    index = find(~cellfun(@isempty,cellSequences)); % non empty cells in cell sequences array
    
    
    
    for j = 1:size(index)
        
        
        
    end
    
    
    clearvars -except files num_files k
end
