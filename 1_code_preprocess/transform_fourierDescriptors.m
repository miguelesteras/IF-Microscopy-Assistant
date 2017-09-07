%% Phase 4.1. Transform cell mask into Fourier descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell mask into Fourier descriptors.
%   Fourier descriptors as first described by Kuhl and Giardina in 
%   "Elliptic Fourier features of a closed contour" 
%   Computer Graphics and Image Processing 18:236-258 1982
%   
%   This code is a modification of: 
%   http://geekstack.net/index.php?page=fourier_descriptor
%   Copyright © Tobias Pohlen 2017
%   A basic implementation of the descriptor presented in:
%   O. D. Trier, A. K. Jain, T. Taxt, Feature Extraction Methods for 
%   Character Recognition - A Survey, Pattern Recognition, Vol. 29, No. 4, 
%   pp. 641 - 662 (1996).
%   ======================================================================

% seqLength = Number of frames in learning sequence (x input + 1 target).
% C determines the number of coefficients (and the description vector
% size). P determines the number of sampling points used for the fourier
% transform (in percentage of total number of points). The transformation
% is set not to be rotation invariant.
seqLength = 8;              
C = 4;                 
P = 0.5;          
rotationInva = false;
%scaleFactor = [linspace(1,2,seqLength-1) 1];
fourierDescriptor = [];
%fourierDescriptorScaled = [];

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    fouriertemp = cell(100,seqLength+1); 
    fouriertempScaled = cell(100,seqLength+1);
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
                % including a change in the coordenates system, from origin
                % [0,0] being center of image to [0,0] being top left image
                % corner and postive y axis = rows (hence 'y*-1' is needed)
                x = round(x)+metadata.maxRadious+1; 
                y = round(y*-1)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), y, x);
                canvas(cellIdx) = true;
                canvas = imfill(canvas,'holes');
                [a,b,c,d,~] = ellipticFourierDescriptor(canvas, C, P, rotationInva);
                fouriertemp{count,m} = [a;b;c;d];
                %fouriertempScaled{count,m} = [a;b;c;d] * scaleFactor(m);
            end
            fouriertemp{count,end} = canvas;
            %fouriertempScaled{count,end} = canvas;
            count = count +1;
        end
    end
    %save(strcat(metadata.name,'_fourierDescriptor.mat'),'fourierDescriptor');
    %save(strcat(metadata.name,'_fourierDescriptorScaled.mat'),'fourierDescriptorScaled');

    fourierDescriptor = [fourierDescriptor;fouriertemp];
    %fourierDescriptorScaled = [fourierDescriptorScaled;fouriertempScaled];
end
save('fourierDescriptorC4.mat','fourierDescriptor');
%save('fourierDescriptorScaled.mat','fourierDescriptorScaled');

%% Image Reconstruction Plots

j = 10;
idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)};            
canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
selection = cellCoordinates{idx(j),idx2(1)};
selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);                
[x,y] = pol2cart(selection(:,1),selection(:,2));
% including a change in the coordenates system, from origin
% [0,0] being center of image to [0,0] being top left image
% corner and postive y axis = rows (hence 'y*-1' is needed)
x = round(x)+metadata.maxRadious+1; 
y = round(y*-1)+metadata.maxRadious+1;
cellIdx = sub2ind(size(canvas), y, x);
canvas(cellIdx) = true;
canvas = imfill(canvas,'holes');

% T = number of sampling points for reconstruction. s = step size.
% ImgSize = size of reconstructed image [rows, columns]. 
T = 300;
s = 0.5;
ImgSize = size(canvas);
C = [1 2 4 6 8 10 30];
dice = zeros(size(C));

for i = 1:numel(C)
    [a,b,c,d,~] = ellipticFourierDescriptor(canvas, C(i), 0.5, rotationInva);
    I = inverseFourierDescriptor(a,b,c,d,T,s,ImgSize);        
%     [B, L] = bwboundaries(I);
    red = zeros(ImgSize);
    red(I==1) = 255;
    RGBimage = double(cat(3, imcomplement(canvas), imcomplement(canvas), imcomplement(canvas)));
    RGBimage(RGBimage==0) = 0.8;
    RGBimage(:,:,1) = RGBimage(:,:,1) + red;
    RGBimage(:,:,2) = RGBimage(:,:,2) - red;
    RGBimage(:,:,3) = RGBimage(:,:,3) - red;
    figure
    imshow(RGBimage);
    % compute Sørensen-Dice Coefficient between original image and reconstructed 
    dice(i) = 2*nnz(canvas&imfill(I,'holes'))/(nnz(imfill(I,'holes')) + nnz(canvas));     
end


%% Plot GFD for KNN

I1 = GFDDescriptors{10,end};
I2 = GFDDescriptors{13,end};
I3 = GFDDescriptors{30,end};
I4 = GFDDescriptors{32,end};
I5 = GFDDescriptors{80,end};
I6 = GFDDescriptors{82,end};

figure;
subplot(2,3,1), imshow(I1)
subplot(2,3,2), imshow(I2)
subplot(2,3,3), imshow(I3)
subplot(2,3,4), imshow(I4)
subplot(2,3,5), imshow(I5)
subplot(2,3,6), imshow(I6)

cell1 = GFDDescriptors{10,1:end-2};
cell2 = GFDDescriptors{13,1:end-2};
cell3 = GFDDescriptors{30,1:end-2};
cell4 = GFDDescriptors{32,1:end-2};
cell5 = GFDDescriptors{80,1:end-2};
cell6 = GFDDescriptors{82,1:end-2};

% Compute Pairwise Euclidean distance between morphologies
dist = pdist([cell1,cell2,cell3,cell4,cell5,cell6]','euclidean');
distSqrF = squareform(dist);

%visualize results
cmap = [0.2422,    0.1504,    0.6603
        0.2803,    0.2782,    0.9221
        0.2440,    0.3358,    0.9988
        0.2120,    0.6314,    0.6901
        0.1709,    0.7154,    0.6342
        0.1938,    0.7758,    0.6251
        0.5044,    0.7993,    0.3480
        0.8634,    0.7406,    0.1596
        0.9892,    0.8136,    0.1885];   
figure;
imagesc(distSqrF);
colormap(cmap);
txt = num2str(distSqrF(:),'%0.2f');
txt = strtrim(cellstr(txt));
[x,y] = meshgrid(1:6);
hstr  = text(x(:),y(:),txt(:),'HorizontalAlignment','center','color','white');
names = {'cell 1','cell 2','cell 3','cell 4','cell 5','cell 6'};
set(gca,'FontSize', 16,'XTickLabel',names,'YTickLabel',names)



