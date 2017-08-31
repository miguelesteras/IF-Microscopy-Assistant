%% Phase 5.7. KNN Search with Generic Fourier descriptors (GFD)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. KNN search (euclidean distance) over a dictionary of GFD of cell sequences.  
%   ======================================================================

%% K-Nearest Neighbours search using stored GFD sequences
rng('default')

load('GFDDescriptors.mat','GFDDescriptors');     
load('GFDDescriptorsScaled.mat','GFDDescriptorsScaled');     

% create dictionary and query set.
queryRatio = .15;
queryIdx = randi([1 size(X,1)],1,round(size(X,1)*queryRatio));
Query = GFDDescriptorsScaled(queryIdx,:); 
Dicti = GFDDescriptorsScaled;
Dicti(queryIdx,:) = [];

X = cell2mat(cellfun(@transpose,Dicti(:,1:end-2),'UniformOutput',0));
Y = cell2mat(cellfun(@transpose,Query(:,1:end-2),'UniformOutput',0));

Ximg = Dicti(:,end);
Yimg = Query(:,end);

% look up amoung all stored sequences (X) the NN of each observation in Y
[IDX,dist] = knnsearch(X,Y,'NSMethod','exhaustive', 'Distance', 'euclidean');

% retrieve predictions
Yhat = Ximg(IDX);

% compare predicted against target shape


% plot shapes
figure
imshowpair(Yimg{4},Yhat{4},'montage');

figure
toGray = 255 * uint8(Yimg{1});
toGray2 = 255 * uint8(Yhat{1});
image = cat(3, toGray, toGray, toGray);
image(:,:,1) = image(:,:,1) + toGray2;
image(:,:,2) = image(:,:,2) - toGray2;
image(:,:,3) = image(:,:,3) - toGray2;
imshow(image)
