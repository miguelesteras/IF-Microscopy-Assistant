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
% Load all tiff files in working directory and store video frames as cell
% array

files = dir('*.tif');      % read tiff videos from folder
num_files = length(files);

for k = 1:num_files
    tiffInfo = imfinfo(files(k).name);  % tiff file inforamtion
    no_frames = numel(tiffInfo);        % number of images in tiff file
    movie = cell(no_frames,1);           % create cell array to contain video frames
    
    for i = 1:no_frames
        I=imread(files(k).name,i);
        movie{i}=I;    
    end
    
    [~,name,~] = fileparts(files(k).name);
    name(~ismember(name,['0':'9' 'a':'z' 'A':'Z'])) = '';  % Remove characters using ismember
    save(strcat(name,'_movie.mat'),'movie');               % save movie to disk as matlab variable 
    
    metadata = struct('name',name,...       % record experiment metadata
                    'imageSize', [],...
                    'minCellArea', [],...    
                    'maxCellArea', [],...
                    'cellMedian', [],...
                    'cellNumbers', [],...
                    'minNucArea', [],...
                    'maxNucArea', [],...
                    'nucleusMedian', []);  
                  
    save(strcat(name,'_metadata.mat'),'metadata');  % save metadata variable 
    
end

clear; clc;

%% Create Binary Masks for Nucleus and Cells

files = dir('*_metadata.mat');      % detect cell arrays created in previous section
num_files = length(files);
for k = 1:num_files
    load(files(k).name,'metadata');    
    load(strcat(metadata.name,'_movie.mat'),'movie');            % load cell array with movie
    
    %% Detect nucleus and create a binary mask

    red_channel = cell(size(movie,1),1);    % cell array containing only the red channel as binary
    nucleus_count = cell(size(movie,1),1);  % record of number of nucleus detected and size
    
    for i = 1:size(movie,1)
        red = movie{i}(:,:,1);      % red channel (where nucleus are)
        red = imbinarize(red,'adaptive',...     % binarize using Otsu's threshold method
            'ForegroundPolarity','bright','Sensitivity',0.1); 
        se = strel('disk',5);
        red = imopen(imfill(red,'holes'), se);     % fill holes and morphological opening
        
        nucleus = bwconncomp(red);              % detect connected components
        red_channel{i} = red;       
        nucleus_count{i} = cellfun('length',nucleus.PixelIdxList);

    end
    
    nucAreas = cell2mat(nucleus_count');
    nucAreas = sort(nucAreas);        % sort nucleus area values 
    idx = round(size(nucAreas,2)/10);
    minNucArea = nucAreas(3*idx);             % set min area of nucleus 
    maxNucArea = nucAreas(end);               % set max area of nucleus (optional)
    nucleusMedian = median(nucAreas(3*idx:end));
    %histogram(nucAreas)
    metadata.imageSize = nucleus.ImageSize; % add info to metadata file
    metadata.minNucArea = minNucArea;    
    metadata.maxNucArea = maxNucArea;
    metadata.nucleusMedian = nucleusMedian;
    
    nucleusInfo = cell(size(red_channel,1),1);      % cell array with nucleus area and location info
    nucleusMask = cell(size(red_channel,1),1);      % cell array with nucleus mask
    nucleusContour = cell(size(red_channel,1),1);   % cell array with nucleus contour
    
    for j = 1:size(red_channel,1)
        frame = bwareafilt(red_channel{j},[minNucArea maxNucArea]);  % remove small and large objects
        contour = bwmorph(frame,'remove');
        nucleusInfo{j} = bwconncomp(frame);         % detect connected components
        nucleusMask{j} = frame;                     % store nulceus binary mask 
        nucleusContour{j} = contour;                % store nucleus binary contour
    end
    
    save(strcat(metadata.name,'_nucleusInfo.mat'),'nucleusInfo');   % save variables in disk
    save(strcat(metadata.name,'_nucleusMask.mat'),'nucleusMask');
    save(strcat(metadata.name,'_nucleusContour.mat'),'nucleusContour');    
    
    %% Detect cells and create a binary mask
       
    cellInfo        = cell(size(movie,1),1);    % initialize record for cells and region properties
    singleCellArea  = cell(size(movie,1),1);    % initialize record for areas of single cells
    cellMask        = cell(size(movie,1),1);    % initialize mask containing all cells in frame
    singleCellMask  = cell(size(movie,1),1);    % initialize mask for single cells
    clusterMask     = cell(size(movie,1),1);    % initialize mask for clusters of cells
    singleCellContour = cell(size(movie,1),1);  % initialize mask for single cells contour
    clusterContour  = cell(size(movie,1),1);    % initialize mask for cluster of cells contour
    
    for m = 1:size(movie,1)
        green = imbinarize(movie{m}(:,:,2),...  % binarize using Otsu's threshold method
            'adaptive',...     
            'ForegroundPolarity','bright',...
            'Sensitivity',0.2); 
        
        se = strel('disk',5);
        green = imopen(imfill(green,'holes'), se);      % fill holes and morphological opening
        [r, c] = find(nucleusMask{m}==1);               % list of coordenates of nucleus pixels
        selection = bwselect(green,c,r,8);              % cellmask regions overlapping with nucleus
        selection = imfill(selection,'holes');          % fill holes in each binary region
        cellMask{m} = selection;                            % record cell mask
        
        labelNucleus = labelmatrix(nucleusInfo{m});             % label nucleus mask
        G = uint8(selection);                                   % convert binary to grayscale
        G(nucleusMask{m}) = (labelNucleus(nucleusMask{m}))+1;   % add nucleus labels (1 not included) to cell image
        %collage = label2rgb(G);                                % collage segmented cells + nucleus
        
        CC = bwconncomp(selection);
        cellInfo{m} = CC;                                   % record cell info
        labelCells = labelmatrix(CC);                       % label cell mask
        
        singleCellMask{m} = false(size(selection)); % intialize masks as empty binary image
        clusterMask{m} = false(size(selection));
        
        for p = 1:max(max(labelCells))              
            if size(unique(G(labelCells==p)),1) == 2            % a single nucleus in connected cell region       
               singleCellMask{m} = or(singleCellMask{m},labelCells==p);
            elseif size(unique(G(labelCells==p)),1) > 2         % more than two nucleus in connected cell region
               clusterMask{m} = or(clusterMask{m},labelCells==p);
            end 
            
        sContour = bwmorph(singleCellMask{m},'remove');     % record contour of cells
        singleCellContour{m} = sContour;
        CCC = bwconncomp(singleCellMask{m});
        singleCellArea{m} = cellfun('length',CCC.PixelIdxList);    % record single cell area info

        cContour = bwmorph(clusterMask{m},'remove');
        clusterContour{m} = cContour;
        end        
    end
    
    save(strcat(metadata.name,'_cellMask.mat'),'cellMask');     % save variables in disk
    save(strcat(metadata.name,'_cellInfo.mat'),'cellInfo');
    save(strcat(metadata.name,'_singleCellArea.mat'),'singleCellArea');
    save(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    save(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');
    save(strcat(metadata.name,'_singleCellContour.mat'),'singleCellContour');
    save(strcat(metadata.name,'_clusterContour.mat'),'clusterContour');
    
    %% average cell size

    cellAreas = cell2mat(singleCellArea');
    cellAreas = sort(cellAreas);         
    idx = round(size(cellAreas,2)/10);
    cellAreas = cellAreas(idx:end-idx);     % remove 10% smallest and biggest
    [N,edges] = histcounts(cellAreas,6);    % put areas in 6 bins
    [~,id] = max(N(4:6));                   % select mode bin
    minCellArea = edges(2+id);              % set min area as the edge of one bin smaller than the mode bin 
    maxCellArea = cellAreas(end);           % set max area of nucleus
    cellMedian = median(cellAreas(sum(N(1:1+id)):end));
    %histogram(cellAreas)
    metadata.minCellArea = minCellArea;
    metadata.maxCellArea = maxCellArea;
    metadata.cellMedian = cellMedian;
    metadata.cellNumbers = sum(N);
                  
    save(strcat(metadata.name,'_metadata.mat'),'metadata');     % save metadata file on disk
    
end

clear; clc;

%% Build data set of single cells

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');

    for j = 1:size(singleCellMask,1)
        
    
    end

    

end
