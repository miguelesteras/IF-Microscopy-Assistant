%% Phase 5.7. Train KNN with Generic Fourier descriptors (GFD)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell contour into GFD. Code by Frederik Kratzert
%   https://uk.mathworks.com/matlabcentral/fileexchange/52643-fd-=-gfd-bw-m-n--implementation-of-the-generic-fourier-descriptors
%   Binary image centered using centerobject funtion by Frederik Kratzert. 
%   ======================================================================

% seqLength = Number of frames in learning sequence (x input + 1 target).
% M is the radial frequency. N is the angular frequency.
seqLength = 8;              
M = 9;                           
N = 10;   
GFDDescriptors = [];
GFDDescriptorsScaled = [];
scaleFactor = [linspace(1,2,seqLength-1) 1];

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    
    % GFdtemp stores descriptor vectors as calculated by the function GFD.
    % In GFDtemp2 the values in the vectors are scaled up in funtion of
    % time with a scale factor in the range [1 2]. 1 for the first input 
    % frame in the series and 2 for the last input frame. The output frame
    % is not scaled up. Hence a KNN classifier using a distance measure
    % like 'euclidean' will weight higher frames closer to the output.
    GFDtemp = cell(100,seqLength+1);
    GFDtemp2 = cell(100,seqLength+1); 

    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};            
            for m = 1:seqLength
                canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
                selection = cellCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);          
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                colSub = round(x)+metadata.maxRadious+1; 
                rowSub = round(y)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                canvas = bwareaopen(imfill(bwmorph(canvas,'bridge'),'holes'),10);                
                canvas = centerobject(canvas);
                descriptors = gfd(canvas,M,N);   
                GFDtemp{count,m} = descriptors(1:end-1);
                GFDtemp2{count,m} = descriptors(1:end-1)*scaleFactor(m);
            end
            % the last column stores the binary image of the target frame
            GFDtemp2{count,end} = canvas;
            GFDtemp{count,end} = canvas;
            count = count +1;
        end
    end
    GFDDescriptors = [GFDDescriptors;GFDtemp];
    GFDDescriptorsScaled = [GFDDescriptorsScaled;GFDtemp2];

    %clearvars -except files num_files i GFDDescriptor seqLength M N
end

save('GFDDescriptors.mat','GFDDescriptors');     
save('GFDDescriptorsScaled.mat','GFDDescriptorsScaled');     

%% K-Nearest Neighbours search using stored GFD sequences

% create training and query set.
X = GFDDescriptorsScaled(:,1:seqLength-1);
X = cell2mat(cellfun(@transpose,X,'UniformOutput',0));
Ximg = GFDDescriptorsScaled(:,end);

testRatio = .15;
testIdx = randi([1 size(X,1)],1,round(size(X,1)*testRatio));
Y = X(testIdx,:);         
Yimg = Ximg(testIdx,:);
X(testIdx,:) = [];
Ximg(testIdx,:) = [];

% look up amoung all stored sequences (X) the NN of each observation in Y
[IDX,D] = knnsearch(X,Y,'NSMethod','exhaustive', 'Distance', 'euclidean');

% retrieve predictions
Yhat = Ximg(IDX);

% compare predicted against target shape


% plot shapes
imshowpair(Yimg{1},Yhat{1},'montage');

toGray = 255 * uint8(Yimg{1});
toGray2 = 255 * uint8(Yhat{1});
image = cat(3, grayImage, grayImage, grayImage);
image(:,:,1) = image(:,:,1) + toGray2;
image(:,:,2) = image(:,:,2) - toGray2;
image(:,:,3) = image(:,:,3) - toGray2;
imshow(image)

