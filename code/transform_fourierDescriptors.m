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

% Number of frames in learning sequence (4 input + 1 target)
seqLength = 5;

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    
    
    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellSequences)),2);
    idx = find(noFrames >= seqLength);
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellSequences(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            frameOne = cellSequences(idx(j),idx(k));
            
            
        end
        
        
    end
    
    

    
    
    index = find(~cellfun(@isempty,cellSequences));
    
           
        
    
    
    clearvars -except files num_files k
end
