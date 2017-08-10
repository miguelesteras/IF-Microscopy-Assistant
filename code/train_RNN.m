%% Phase 5.1 Train a RNN
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a Recurrent Neural Network based on Fourier descriptors. 
%   ======================================================================

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                               
    load(strcat(metadata.name,'_fourierInput.mat'),'fourierInput');     
    load(strcat(metadata.name,'_fourierTarget.mat'),'fourierTarget');
    
end
