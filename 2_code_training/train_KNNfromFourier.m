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
 
% K-Nearest Neighbours search using stored fourier sequences

load(strcat(metadata.name,'fourierDescriptor.mat'),'fourierDescriptor');
load(strcat(metadata.name,'fourierDescriptorScaled.mat'),'fourierDescriptorScaled');

% create training and query set.
X = fourierDescriptorScaled(:,1:seqLength-1);
X = cell2mat(cellfun(@transpose,X,'UniformOutput',0));
Ximg = fourierDescriptorScaled(:,end);

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

% quantify difference between predicted and target shape


% plot shapes
figure
imshowpair(Yimg{2},Yhat{2},'montage');

figure
toGray = 255 * uint8(Yimg{1});
toGray2 = 255 * uint8(Yhat{1});
image = cat(3, toGray, toGray, toGray);
image(:,:,1) = image(:,:,1) + toGray2;
image(:,:,2) = image(:,:,2) - toGray2;
image(:,:,3) = image(:,:,3) - toGray2;
imshow(image)
