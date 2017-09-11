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

% Number of frames in learning sequence (seqLength). Construct square
% grayscale image of size doble the maximum radius of cells (maxRad).
seqLength = 8; 
intScale = linspace(20,120,seqLength-1);
dynamicImages = [];
maxRad = 40;

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    tempDynamic = cell(100,2); 
    count = 1;
    
    % the black canvas where to draw the dynamic image has dimensions given
    % by the max radious detected in sample, +1 pixel. This makes is square
    % and containing a central pixel (center of mass).
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};
            input = int8(zeros((maxRad*2)+1,(maxRad*2)+1));
            m = 1;            
            while m < seqLength                
                canvas = false(size(input));
                selection = cellCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);          
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                % including a change in the coordenates system, from origin
                % [0,0] being center of image to [0,0] being top left image
                % corner and postive y axis = rows (hence 'y*-1' is needed)
                colSub = round(x)+(maxRad+1); 
                rowSub = round(y*-1)+(maxRad+1);
                colSub(colSub>(maxRad*2)+1) = (maxRad*2);
                rowSub(rowSub>(maxRad*2)+1) = (maxRad*2);
                colSub(colSub<1) = 1;
                rowSub(rowSub<1) = 1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                SE = strel('disk',2);
                contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
                input(contour) = intScale(m);
                m = m+1;
            end
            tempDynamic{count,1} = mat2gray(input);
            target = int8(zeros(size(input)));
            canvas = false(size(input));
            selection = cellCoordinates{idx(j),idx2(k+m-1)};
            selection(:,1) = selection(:,1) + rotation;          
            [x,y] = pol2cart(selection(:,1),selection(:,2));
            colSub = round(x)+(maxRad+1); 
            rowSub = round(y*-1)+(maxRad+1);
            colSub(colSub>(maxRad*2)+1) = (maxRad*2)+1;
            rowSub(rowSub>(maxRad*2)+1) = (maxRad*2)+1;
            colSub(colSub<1) = 1;
            rowSub(rowSub<1) = 1;
            cellIdx = sub2ind(size(canvas), rowSub, colSub);
            canvas(cellIdx) = true;
            SE = strel('disk',2);
            contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
            target(contour) = 255;
            tempDynamic{count,2} = mat2gray(target);
            count = count+1;
        end               
    end
    dynamicImages = [dynamicImages ; tempDynamic];
    save(strcat(metadata.name,'_dynamicImages.mat'),'tempDynamic');
end
save('dynamicImages.mat','dynamicImages');     


%% cell areas for data augmentation

sequenceArea = [];
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    tempArea = zeros(100,1); 
    count = 1;
    
    % the black canvas where to draw the dynamic image has dimensions given
    % by the max radious detected in sample, +1 pixel. This makes is square
    % and containing a central pixel (center of mass).
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            area = size(cellCoordinates{idx(j),idx2(k)},1);
            tempArea(count) = area;
            count = count+1;
        end
    end
    sequenceArea = [sequenceArea ; tempArea];
end
save('sequenceArea.mat','sequenceArea');     

quant = quantile(sequenceArea,[0.15 0.35 0.65 0.9]);
figure;
histfit(sequenceArea); hold on

Xline = quant;
for idx = 1 : numel(Xline)
    plot([Xline(idx) Xline(idx)], [0 130]);
end


%% Plots

intScale = linspace(20,140,seqLength-1);
input = int8(zeros((maxRad*2)+1,(maxRad*2)+1));
m = 1; 
k = 1;
while m < seqLength                
    canvas = false(size(input));
    selection = cellCoordinates{idx(j),idx2(k+m-1)};
    selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);          
    [x,y] = pol2cart(selection(:,1),selection(:,2));
    % including a change in the coordenates system, from origin
    % [0,0] being center of image to [0,0] being top left image
    % corner and postive y axis = rows (hence 'y*-1' is needed)
    colSub = round(x)+(maxRad+1); 
    rowSub = round(y*-1)+(maxRad+1);
    colSub(colSub>(maxRad*2)+1) = (maxRad*2);
    rowSub(rowSub>(maxRad*2)+1) = (maxRad*2);
    colSub(colSub<1) = 1;
    rowSub(rowSub<1) = 1;
    cellIdx = sub2ind(size(canvas), rowSub, colSub);
    canvas(cellIdx) = true;
    SE = strel('disk',2);
    contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
    figure;
    imshow(imresize(contour,4));
    input(contour) = intScale(m);
    I = mat2gray(input);
    rgb = cat(3, I, I, I);
    figure;
    imshow(imresize(rgb,4))
    m = m+2;
end


%% Plots

No = [200 20 650 50 345];
for i = 1:numel(No)
    n = No(i);
    canvas = cat(3, dynamicImages{n,1}, dynamicImages{n,1}, dynamicImages{n,1});
    figure;
    imshow(imresize(canvas,2))
    canvas(:,:,1) = canvas(:,:,1) + dynamicImages{n,2};
    canvas(:,:,2) = canvas(:,:,2) - dynamicImages{n,2};
    canvas(:,:,3) = canvas(:,:,3) - dynamicImages{n,2};
    figure;
    imshow(imresize(canvas,2))
    output = cat(3, dynamicImages{n,2}, zeros(size(dynamicImages{n,2}),'like',dynamicImages{n,2}),...
                                        zeros(size(dynamicImages{n,2}),'like',dynamicImages{n,2}));
    figure;
    imshow(imresize(output,2))
end