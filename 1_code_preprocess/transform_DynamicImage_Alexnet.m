%% Phase 4.2.1 Transform cell silhouettes into Dynamic images for Alexnet
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
intScale = linspace(20,200,seqLength-1);

% columns 1-2 are original data; other columns is synthetic 
dynamicImagesT = [];

load('metadata.mat','metadata');
maxRad = 113;

% Calculate area values for cummulative quantiles, later used for training
% data synthesis
load('singleCellArea.mat','singleCellArea');
cellAreas = sort(cell2mat(singleCellArea'));
QQ = quantile(cellAreas,[0.4 0.5 0.6 0.7]);

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    
    tempDynamic = cell(100,2*(numel(QQ)+1));    
    count = 1;
    % the black canvas where to draw the dynamic image has dimensions given
    % by the max radious detected in sample, +1 pixel. This makes is square
    % and containing a central pixel (center of mass).
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};
            SeqArea = size(cellCoordinates{idx(j),idx2(k)},1);
            FF = [1 sqrt(QQ/SeqArea)];
            for s = 1:numel(FF)
                input = int8(zeros((maxRad*2)+1,(maxRad*2)+1));                
                m = 1;            
                while m < seqLength    
                    canvas = false(size(input));
                    selection = cellCoordinates{idx(j),idx2(k+m-1)};
                    selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
                    selection(:,2) = selection(:,2)*FF(s);

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
                tempDynamic{count,(s*2)-1} = mat2gray(input);
                tempDynamic{count,s*2} = mat2gray(target);
            end
            count = count+1;
        end               
    end
    dynamicImagesT = [dynamicImagesT ; tempDynamic];
end

dynamicImages  = dynamicImagesT(:,1);
dynamicImages04 = dynamicImagesT(:,3);
dynamicImages05 = dynamicImagesT(:,5);
dynamicImages06 = dynamicImagesT(:,7);
dynamicImages07 = dynamicImagesT(:,9);

save('dynamicImages_Alexnet.mat','dynamicImages');     
save('dynamicImages_Alexnet04.mat','dynamicImages04');     
save('dynamicImages_Alexnet05.mat','dynamicImages05');     
save('dynamicImages_Alexnet06.mat','dynamicImages06');     
save('dynamicImages_Alexnet07.mat','dynamicImages07');     
