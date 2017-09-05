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

% Number of frames in learning sequence (x input + 1 target)
seqLength = 8; 
dynamicImages = [];

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
            input = int8(zeros((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1));
            m = 1;            
            while m < seqLength                
                canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
                selection = cellCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);          
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                % including a change in the coordenates system, from origin
                % [0,0] being center of image to [0,0] being top left image
                % corner and postive y axis = rows (hence 'y*-1' is needed)
                colSub = round(x)+metadata.maxRadious+1; 
                rowSub = round(y*-1)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                SE = strel('disk',2);
                contour = bwmorph(imopen(imfill(canvas,'holes'),SE),'remove');
                input(contour) = 15*m;
                m = m+1;
            end
            tempDynamic{count,1} = mat2gray(input);
            target = int8(zeros((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1));
            canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
            selection = cellCoordinates{idx(j),idx2(k+m-1)};
            selection(:,1) = selection(:,1) + rotation;          
            [x,y] = pol2cart(selection(:,1),selection(:,2));
            colSub = round(x)+metadata.maxRadious+1; 
            rowSub = round(y*-1)+metadata.maxRadious+1;
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

%% Plots

canvas = cat(3, dynamicImages{5,1}, dynamicImages{5,1}, dynamicImages{5,1});
canvas(:,:,1) = canvas(:,:,1) + dynamicImages{5,2};
canvas(:,:,2) = canvas(:,:,2) - dynamicImages{5,2};
canvas(:,:,3) = canvas(:,:,3) - dynamicImages{5,2};
imshow(canvas)

%% Transform into gradients

% [Gmag,Gdir] = imgradient(dynamicInput{1});
% [Gx, Gy] = imgradientxy(dynamicInput{1},'prewitt');
% figure
% imshowpair(Gx, Gy,'montage');
% C = imfuse(Gx, Gy, 'montage');
% figure
% imshowpair(Gmag,Gdir,'montage');