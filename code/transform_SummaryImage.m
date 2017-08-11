%% Phase 4.2. Transform cell silhouettes into Dynamic images
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell silhouettes of input sequence into a single gray scale 
%   image (summary of all input images).
%
%   ======================================================================

%% Build data set of single cells

seqLength = 5;              % Number of frames in learning sequence (4 input + 1 target)

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_centerMass.mat'),'centerMass');     
    
    cellCoordenates = cell(size(cellSequences));
    index = find(~cellfun(@isempty,cellSequences)); % non empty cells in cell sequences array
    for h = 1:size(index)    
        selection = cellSequences{index(h)};
        center = centerMass{index(h)};
        image = false(metadata.imageSize);
        image(selection) = true;                % only show selected cell in binary image
        stats = regionprops(image,'PixelList');
        coordenates = stats.PixelList - center;
        [theta, rho] = cart2pol(coordenates(:,1),coordenates(:,2));   % polar coordenates
        cellCoordenates{index(h)} = [theta rho]; 
    end

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordenates)),2);
    idx = find(noFrames >= seqLength);
    dynamicInput = cell(100,1); 
    dynamicTarget = cell(100,1);
    count = 1;
    
    % the black canvas where to draw the dynamic image has dimensions given
    % by the max radious detected in sample, +1 pixel. This makes is square
    % and containing a central pixel (center of mass).
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordenates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};
            input = int8(zeros((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1));
            m = 1;            
            while m < seqLength                
                canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
                selection = cellCoordenates{idx(j),idx2(k+m-1)};
                selection(:,1) = selection(:,1) + rotation;          % apply rotation
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                colSub = round(x)+metadata.maxRadious+1; 
                rowSub = round(y)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                SE = strel('disk',2);
                contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
                input(contour) = 50*m;
                m = m+1;
            end
            dynamicInput{count,1} = mat2gray(input);
            target = int8(zeros((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1));
            canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
            selection = cellCoordenates{idx(j),idx2(k+m-1)};
            selection(:,1) = selection(:,1) + rotation;          
            [x,y] = pol2cart(selection(:,1),selection(:,2));
            colSub = round(x)+metadata.maxRadious+1; 
            rowSub = round(y)+metadata.maxRadious+1;
            cellIdx = sub2ind(size(canvas), rowSub, colSub);
            canvas(cellIdx) = true;
            SE = strel('disk',2);
            contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
            target(contour) = 50*m;
            dynamicTarget{count,1} = mat2gray(target);
            count = count+1;
        end               
    end
    save(strcat(metadata.name,'_dynamicInput.mat'),'dynamicInput');     
    save(strcat(metadata.name,'_dynamicTarget.mat'),'dynamicTarget');     
    save(strcat(metadata.name,'_cellCoordenates.mat'),'cellCoordenates');     
    
    clearvars -except files num_files i
end



load(files(1).name,'metadata');                               
load(strcat(metadata.name,'_fourierInput.mat'),'fourierInput');     
load(strcat(metadata.name,'_fourierTarget.mat'),'fourierTarget');
timeSeries = [fourierInput fourierTarget];   

for i = 2:num_files
    load(files(i).name,'metadata');                               
    load(strcat(metadata.name,'_fourierInput.mat'),'fourierInput');     
    load(strcat(metadata.name,'_fourierTarget.mat'),'fourierTarget');
    
    timeSeries = [timeSeries; [fourierInput fourierTarget]];
end
