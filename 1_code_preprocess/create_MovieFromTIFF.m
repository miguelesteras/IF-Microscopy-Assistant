%% Phase 1.0 Transform TIFF to cell array Create Masks
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================%% Pre-processing;
% All .tiff files must be placed in working directory.
% They will be transformed into individual cell arrays containing the video
% frames

files = dir('*.tif');      
num_files = length(files);

for k = 1:num_files
    % extract tiff file information, including number of frames. Create an
    % empty cell array to store all video frames
    tiffInfo = imfinfo(files(k).name);  
    no_frames = numel(tiffInfo);        
    movie = cell(no_frames,1);           
    
    for i = 1:no_frames
        I=imread(files(k).name,i);
        movie{i}=I;    
    end
    
    % save movie cell array to disk as matlab variable.  
    [~,name,~] = fileparts(files(k).name);
    name(~ismember(name,['0':'9' 'a':'z' 'A':'Z'])) = '';
    save(strcat(name,'_movie.mat'),'movie');               
    
    % record experiment metadata and save it as matlab variable
    metadata = struct('name',name,...       
                    'imageSize', [],...
                    'minCellArea', [],...    
                    'maxCellArea', [],...
                    'cellMedian', [],...
                    'cellNumbers', [],...
                    'minNucArea', [],...
                    'maxNucArea', [],...
                    'nucleusMedian', []);                
    save(strcat(name,'_metadata.mat'),'metadata');
    
    clearvars -except files num_files k
end
clear; clc;
