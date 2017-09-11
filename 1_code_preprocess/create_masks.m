%% Phase 1. Create Masks
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Detect nucleus           -> create binary mask for nucleus
%   2. Detect single cells      -> create binary mask for cells
%   3. Detect clusters of cells -> create binary mask for clusters
%   ======================================================================

%% Pre-processing;
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

%% Reduce Image Size (Gaussian pyramid) and Create Binary Masks for Nucleus and Cells

files = dir('*_metadata.mat');
num_files = length(files);
for k = 1:num_files
    load(files(k).name,'metadata');    
    load(strcat(metadata.name,'_movie.mat'),'movie'); 
    
    %% Compute a Gaussian pyramid reduction by one level
    for j = 1:size(movie,1)
        movie{j} = impyramid(movie{j}, 'reduce');    
    end

    %% Detect nucleus and create a binary mask
    % the red_channel array will store the binary mask of the red channel
    % where the nucleus are detected. The array nucleus_count will
    % record the number of nucleus and size found in each frame.
    red_channel = cell(size(movie,1),1);    
    nucleus_count = cell(size(movie,1),1);
    
    for i = 1:size(movie,1)
        red = movie{i}(:,:,1);
        % binarize the red channel image using Otsu's threshold method
        red = imbinarize(red,'adaptive','ForegroundPolarity','bright','Sensitivity',0.1); 
        % fill holes and morphological opening
        se = strel('disk',5);
        red = imopen(imfill(red,'holes'), se);     
        % Properties of the detect connected components are stored
        nucleus = bwconncomp(red);              
        red_channel{i} = red;       
        nucleus_count{i} = cellfun('length',nucleus.PixelIdxList);
    end
    % set min and max area of nucleus (20% smallest removed) and store
    % information on metadata file. 
    nucAreas = sort(cell2mat(nucleus_count')); 
    idx = round(size(nucAreas,2)/10);
    minNucArea = nucAreas(2*idx);              
    maxNucArea = nucAreas(end);               
    nucleusMedian = median(nucAreas(3*idx:end));
    %histogram(nucAreas)
    metadata.imageSize = nucleus.ImageSize;
    metadata.minNucArea = minNucArea;    
    metadata.maxNucArea = maxNucArea;
    metadata.nucleusMedian = nucleusMedian;
    
    % create cell arrays with nucleus count, area and location
    % (nucleusInfo), and nucleus binary mask (nucleusMask).
    nucleusInfo = cell(size(red_channel,1),1);     
    nucleusMask = cell(size(red_channel,1),1);      
    
    % nucleus included in binary mask are thresholded by min and max values
    % calculated before.
    for j = 1:size(red_channel,1)
        frame = bwareafilt(red_channel{j},[minNucArea maxNucArea]);  
        nucleusInfo{j} = bwconncomp(frame);         
        nucleusMask{j} = frame;                      
    end
    
    %% Detect cells and create a binary mask
    
    % Initialize cell arrays of binary mask for all cell
    % structures (cellMask), only single cells (singleCellMask), only cell
    % clusters (clusterMask), and cell arrays to contain single cell areas
    % (singleCellAreas). 
    singleCellArea  = cell(size(movie,1),1);   
    cellMask        = cell(size(movie,1),1);   
    singleCellMask  = cell(size(movie,1),1);   
    clusterMask     = cell(size(movie,1),1);   

    for m = 1:size(movie,1)
        % hysteresis thresholding. Threshold level found using using 
        % Otsu's method
        [level,~] = graythresh(movie{m}(:,:,2));
        [green,~] = hysteresis3d(movie{m}(:,:,2),level*0.7,level,4);        
        % fill holes and morphological opening
        se = strel('disk',4);
        green = imopen(imfill(green,'holes'), se);
        % foreground regions in the binary cellmask regions that overlap
        % with foreground pixels in the binary nucleus mask are selected as
        % containing a cell structure
        [r, c] = find(nucleusMask{m}==1);               
        selection = bwselect(green,c,r,8);              
        selection = imfill(selection,'holes');          
        cellMask{m} = selection;  
        CC = bwconncomp(selection);
        labelCells = labelmatrix(CC);                       
        
        % label nucleus mask, and add labels to cell mask to create a
        % collage of segmented cells and nucleus. This composed image will
        % be used to divide single cells from cell clusters
        labelNucleus = labelmatrix(nucleusInfo{m});             
        G = uint8(selection);                                   
        G(nucleusMask{m}) = (labelNucleus(nucleusMask{m}))+1;   
        %collage = label2rgb(G);                             
              
        % intialize masks as empty binary image
        singleCellMask{m} = false(size(selection)); 
        clusterMask{m} = false(size(selection));
        % a single nucleus in a connected cell region = single cell
        % more than two nucleus in connected cell region = cell cluster
        for p = 1:max(max(labelCells))              
            if size(unique(G(labelCells==p)),1) == 2                 
               singleCellMask{m} = or(singleCellMask{m},labelCells==p);
            elseif size(unique(G(labelCells==p)),1) > 2         
               clusterMask{m} = or(clusterMask{m},labelCells==p);
            end       
        % record single cell area info. Later used to identified cell in
        % partial or full view (due to changes in depth). 
        CCC = bwconncomp(singleCellMask{m});
        singleCellArea{m} = cellfun('length',CCC.PixelIdxList);   
        end        
    end
    % save variables into disk
    save(strcat(metadata.name,'_movieReduced.mat'),'movie');
    save(strcat(metadata.name,'_nucleusInfo.mat'),'nucleusInfo');
    save(strcat(metadata.name,'_nucleusMask.mat'),'nucleusMask');
    save(strcat(metadata.name,'_cellMask.mat'),'cellMask');     
    save(strcat(metadata.name,'_singleCellArea.mat'),'singleCellArea');
    save(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    save(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');
    
    %% average cell size

    % remove outliers (10% smallest and biggest), then select mode bin from
    % biggest half. Set min area as the edge of two bins smaller from the
    % mode bin.
    cellAreas = sort(cell2mat(singleCellArea'));         
    idx = round(size(cellAreas,2)/10);
    cellAreas = cellAreas(idx:end-idx);     
    [N,edges] = histcounts(cellAreas,8);    
    [~,id] = max(N(5:8));                   
    minCellArea = edges(2+id);               
    maxCellArea = cellAreas(end);           
    cellMedian = median(cellAreas(sum(N(1:1+id))+1:end));
    %histogram(cellAreas)
    metadata.minCellArea = minCellArea;
    metadata.maxCellArea = maxCellArea;
    metadata.cellMedian = cellMedian;
    metadata.cellNumbers = sum(N);              
    save(strcat(metadata.name,'_metadata.mat'),'metadata'); 
    
    clearvars -except files num_files k
end
