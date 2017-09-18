%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Load variables

files = dir('*_metadata.mat');
load(files(2).name,'metadata');    
load(strcat(metadata.name,'_movie.mat'),'movie'); 
    

%% Plots nucleus

frame = movie{30};
imshow(frame*1.8)
I = frame(:,:,1);
[level,~] = graythresh(I);

figure
imhist(I); hold on;
Xline = level*255;
for idx = 1 : numel(Xline)
    plot([Xline(idx) Xline(idx)], [0 250]);
end
hold off;

figure
imshow(cat(3, I*2, ones(size(I)), ones(size(I))));

IB = imbinarize(I);
se = strel('disk',2);
IB = imopen(imfill(IB,'holes'), se);
figure
imshow(IB)

IB2 = imbinarize(I,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
se = strel('disk',2);
IB2 = imopen(imfill(IB2,'holes'), se);
figure;
imshow(IB2)

[IB3,~] = hysteresis3d(I,level*0.5,level,4);
se = strel('disk',2);
IB3 = imopen(imfill(IB3,'holes'), se);
figure
imshow(IB3)


%% Plots single vs hysteresis thresholding for cells

frame = movie{80};
imshow(frame*1.8)
I = frame(:,:,2);
[level,~] = graythresh(I);

figure
imshow(cat(3, ones(size(I)), I*2.2, ones(size(I))));

IB = imbinarize(I);
se = strel('disk',2);
IB = imopen(imfill(IB,'holes'), se);
figure
imshow(IB)

IB2 = imbinarize(I,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4); 
se = strel('disk',2);
IB2 = imopen(imfill(IB2,'holes'), se);
figure
imshow(IB2)

[IB3,hys] = hysteresis3d(I,level*0.5,level*0.9,4);
se = strel('disk',3);
IB3 = imclose(imopen(imfill(IB3,'holes'), se),se);
figure
imshow(IB3)

%% Plots on morphological opening and imfill

I = movie{50};
figure; imshow(I*1.8)

green = I(:,:,2);
figure; imshow(cat(3, ones(size(green)), green*2.2, ones(size(green))));

[level,~] = graythresh(green);
[BW,~] = hysteresis3d(green,level*0.7,level,4);
figure; imshow(BW)

BW2 = imfill(BW,'holes');
figure; imshow(BW2)

se = strel('disk',4);
BW3 = imopen(BW2, se);
figure; imshow(BW3)

%% Plot on single vs cluster cells

I = movie{400};
imshow(I*1.8)
red = I(:,:,1);
green = I(:,:,2);

figure; imshow(cat(3, ones(size(green)), green*2.2, ones(size(green))));
figure; imshow(cat(3, red*2, ones(size(red)), ones(size(red))));

IB = imbinarize(red,'adaptive','ForegroundPolarity','bright','Sensitivity',0.2);
se = strel('disk',2);
redBW = imclose(imopen(imfill(IB,'holes'), se),se);
figure; imshow(redBW)

[level,~] = graythresh(green);
[BW,~] = hysteresis3d(green,level*0.7,level,4);
se = strel('disk',4);
greenBW = imopen(imfill(BW,'holes'), se);
figure; imshow(greenBW)

[r, c] = find(redBW==1);               
selection = bwselect(greenBW,c,r,8);              
selection = imfill(selection,'holes');          
CC = bwconncomp(selection);
labelCells = labelmatrix(CC);                       

% label nucleus mask, and add labels to cell mask to create a
% collage of segmented cells and nucleus. This composed image will
% be used to divide single cells from cell clusters
nucleusInfo = bwconncomp(redBW);         
labelNucleus = labelmatrix(nucleusInfo); 

G = uint8(selection);    
G(selection) = 2;
G(redBW) = (labelNucleus(redBW))+6;   
collage = label2rgb(G);
figure; imshow(collage);

singleCellMask = false(size(selection)); 
clusterMask = false(size(selection));
for p = 1:max(max(labelCells))              
    if size(unique(G(labelCells==p)),1) == 2                 
       singleCellMask = or(singleCellMask,labelCells==p);
    elseif size(unique(G(labelCells==p)),1) > 2         
       clusterMask = or(clusterMask,labelCells==p);
    end       
end        
figure; imshow(singleCellMask);
figure; imshow(clusterMask);

%% Plots for cell areas

load('singleCellArea.mat','singleCellArea');
cellAreas = sort(cell2mat(singleCellArea'));
count = floor(numel(cellAreas)/10);
figure; histogram(cellAreas,count, 'FaceColor',[.7 .7 .7], 'EdgeColor', [.85 .85 .85]); 
limits = quantile(cellAreas,[0.001 0.99]);  
xlim([limits(1) limits(2)])
xlabel('Single Cell Area (pixels)')
ylabel('Count')
set(gca,'fontsize',16)
hold on;
lines = quantile(cellAreas,[0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9]);  

yl = ylim;
plot([median(cellAreas) median(cellAreas)], [0 yl(2)], 'k', 'LineWidth',2);

for idx = 1 : numel(lines)
    plot([lines(idx) lines(idx)], [0 yl(2)], 'r', 'LineWidth',1);
end
hold off;

% idx = round(size(cellAreas,2)/10);
% cellAreas2 = cellAreas(idx:end-idx);     
% [N,edges] = histcounts(cellAreas2,8);
% [~,id] = max(N(5:8));                   
% minCellArea = edges(2+id);               
% maxCellArea = cellAreas2(end);           
% cellMedian = median(cellAreas2(sum(N(1:1+id))+1:end));
% figure; histogram(cellAreas2,8)

%% Filter single cell area

files = dir('*_metadata.mat');      
OneMetadata = load(files(3).name,'metadata');
OneMetadata = OneMetadata.metadata;
load(strcat(OneMetadata.name,'_singleCellMask.mat'),'singleCellMask');
load(strcat(OneMetadata.name,'_movieReduced.mat'),'movie');       
load(strcat(OneMetadata.name,'_cellMask.mat'),'cellMask');       
load('metadata.mat','metadata');      


j = 30;

I = movie{j};
imshow(I*1.8)
red = I(:,:,1);
green = I(:,:,2);

figure; imshow(cat(3, ones(size(green)), green*2.2, ones(size(green))));
figure; imshow(cat(3, red*2, ones(size(red)), ones(size(red))));

IB = imbinarize(red,'adaptive','ForegroundPolarity','bright','Sensitivity',0.2);
se = strel('disk',2);
redBW = imclose(imopen(imfill(IB,'holes'), se),se);
figure; imshow(redBW)

[level,~] = graythresh(green);
[BW,~] = hysteresis3d(green,level*0.7,level,4);
se = strel('disk',4);
greenBW = imopen(imfill(BW,'holes'), se);
figure; imshow(greenBW)

[r, c] = find(redBW==1);               
selection = bwselect(greenBW,c,r,8);              
selection = imfill(selection,'holes');          
CC = bwconncomp(selection);
labelCells = labelmatrix(CC);                       

nucleusInfo = bwconncomp(redBW);         
labelNucleus = labelmatrix(nucleusInfo); 

G = uint8(selection);    
G(selection) = 2;
G(redBW) = (labelNucleus(redBW))+6;   
collage = label2rgb(G);
figure; imshow(collage);

singleCellMask = false(size(selection)); 
clusterMask = false(size(selection));
for p = 1:max(max(labelCells))              
    if size(unique(G(labelCells==p)),1) == 2                 
       singleCellMask = or(singleCellMask,labelCells==p);
    elseif size(unique(G(labelCells==p)),1) > 2         
       clusterMask = or(clusterMask,labelCells==p);
    end       
end        
figure; imshow(singleCellMask);
figure; imshow(clusterMask);

selection = bwareafilt(singleCellMask,[metadata.minCellArea metadata.maxCellArea]);
figure; imshow(selection)


%% Plots on cell sequences

summaryLength = [];
files = dir('*_metadata.mat');      
load(files(2).name,'metadata');
load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');

% create cell array with the location of only single cell in full
% view as defined to be between min and max cell area.
fullCellLocation = cell(size(singleCellMask,1),1);   

for j = 1:size(singleCellMask,1)
    selection = bwareafilt(singleCellMask{j},[metadata.minCellArea metadata.maxCellArea]); % select only cells in full view  
    pixels = struct2cell(regionprops(selection,'PixelIdxList'));
    fullCellLocation{j} = pixels;     
end


MeanLength = [];
NumSequence = [];
NumSeq8     = [];

X = [0.5 0.6 0.7 0.75 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 0.99];
for i = 1:numel(X)  
    S = X(i);

    % cell array to store cell sequences
    cellSequences = cell(1,size(singleCellMask,1)); 
    seqCount = 1;

    % detect cells and start sequence in first frame
    for m = 1:size(fullCellLocation{1},2)   
        object = fullCellLocation{1}{m};
        cellSequences{seqCount,1} = object;
        seqCount = seqCount + 1;
    end

    % from second frame onwards, for every cell in frame, check if the cell
    % was present in previoius frame (based on >= similar location). If
    % true, add cell to previous sequence, if false, start a new sequence.
    for k = 2:size(fullCellLocation,1)  
        if isempty(fullCellLocation{k})
            continue
        else
            preSeq = find(~cellfun(@isempty,cellSequences(:,k-1)));
            prior = cellSequences(preSeq,k-1);
            for m = 1:size(fullCellLocation{k},2)
                object = fullCellLocation{k}(m);
                object = repmat(object,[size(prior,1) 1]);
                simil = cellfun('length',cellfun(@intersect, prior, object, 'UniformOutput', false));
                [value, idx] = max(simil);                          
                if isempty(simil) || value < size(fullCellLocation{k}{m},1)*S
                    cellSequences{seqCount,k} = fullCellLocation{k}{m};
                    seqCount = seqCount + 1;
                else    
                    cellSequences{preSeq(idx),k} = fullCellLocation{k}{m};                      
                end
            end
        end
    end

    % Binary map of sequences, histogram and sequence lengths
    binarymap = ~cellfun(@isempty,cellSequences);
    %figure; imshow(binarymap);                      
    lengths = sum(uint8(binarymap),2);
    summaryLength = [summaryLength;lengths(lengths > 1)];
    %figure; histogram(summaryLength,90)

    % total number of unique sequences (>1 frames)
    numSeq = size(lengths(lengths > 1),1);
    % total number of 8 frames sequences (potential training data size)
    trainSize = sum(lengths(lengths > 7)-7);  
    
    MeanLength = [MeanLength ; mean(lengths)];
    NumSequence = [NumSequence ; numSeq];
    NumSeq8     = [NumSeq8 ; trainSize];
end

ticks = [0.5 0.6 0.7 0.8 0.9 1];

figure; plot(X,MeanLength, 'Color', 'green', 'LineWidth', 4); hold on
yl = ylim;
plot([0.8 0.8], [0 yl(2)], 'r', 'LineWidth',2);
grid on
xticks(ticks)
xlabel('Similarity')
ylabel('Mean Length of Sequences')
set(gca,'fontsize',16)

figure; plot(X,NumSequence, 'Color', 'green', 'LineWidth', 4); hold on
yl = ylim;
plot([0.8 0.8], [0 yl(2)], 'r', 'LineWidth',2);
grid on
xticks(ticks)
xlabel('Similarity')
ylabel('Total No. of Unique Sequences (>=2 frame)')
set(gca,'fontsize',16)

figure; plot(X,NumSeq8, 'Color', 'green', 'LineWidth', 4); hold on
yl = ylim;
plot([0.8 0.8], [0 yl(2)], 'r', 'LineWidth',2);
grid on
xticks(ticks)
xlabel('Similarity')
ylabel('Total No. of Sequences (>=8 frames)')
set(gca,'fontsize',16)

%% sequence length histogram

len = sort(summaryLength);
quant8 = find(len==8, 1, 'first')/numel(summaryLength);
quant4 = find(len==4, 1, 'first')/numel(summaryLength);

bins = linspace(1,51,51);
figure; histogram(summaryLength,bins,'FaceColor',[.7 .7 .7], 'EdgeColor', [.85 .85 .85]); hold on
yl = ylim;
plot([8 8], [0 yl(2)], 'r', 'LineWidth',2);
xlim([0 50])
ylabel('Count')
xlabel('Unique Sequence Length (frames)')
set(gca,'fontsize',16)

%% rotation normalization

ps1 = polarscatter(wrapTo2Pi(theta),rho, 'filled');
ps1.SizeData = 50;
ps1.MarkerFaceAlpha = .5;
hold on
ps2 = polarscatter(wrapTo2Pi(theta+rotation),rho,'filled');
ps2.SizeData = 50;
ps2.MarkerFaceAlpha = .7;
hold off
lg = legend('Original','Normalized Rotation');
lg.FontSize = 14;

%% Rotation normalization II

files = dir('*_metadata.mat');      
load(files(1).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
idx = find(noFrames >= 8);

j = 2;
idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)};
sequence = cell(1,8);
for m = 1:8
    canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
    selection = cellCoordinates{idx(j),idx2(1+m-1)};
    [x,y] = pol2cart(selection(:,1),selection(:,2));
    % including a change in the coordenates system, from origin
    % [0,0] being center of image to [0,0] being top left image
    % corner and postive y axis = rows (hence 'y*-1' is needed)
    x = round(x)+metadata.maxRadious+1; 
    y = round(y*-1)+metadata.maxRadious+1;
    cellIdx = sub2ind(size(canvas), y, x);
    canvas(cellIdx) = true;
    canvas = imfill(canvas,'holes');
    sequence{m} = canvas;
end

sequenceR = cell(1,8);
for m = 1:8
    canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
    selection = cellCoordinates{idx(j),idx2(1+m-1)};
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
    sequenceR{m} = canvas;
end

for img = 1:numel(sequence)
    subplot(2,8,img); imshow(sequence{img})
    subplot(2,8,8+img); imshow(sequenceR{img})
end

%% Image Reconstruction Plots

files = dir('*_metadata.mat');      
load(files(1).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
idx = find(noFrames >= 8);
        
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


%% Plot GFD VS Elliptic Fourier for KNN

load('GFDDescriptors.mat','GFDDescriptors');
load('FourierDescriptorC10.mat','FourierDescriptorC10');

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

%% Plot rotation in boundary descriptors

files = dir('*_metadata.mat');      
load(files(3).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
idx = find(noFrames >= 8);
        
j = 2;
idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)} + pi/2;            
canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
selection = cellCoordinates{idx(j),idx2(1)};
selection(:,1) = wrapTo2Pi(selection(:,1) + rotation); 
figure;
polarscatter(selection(:,1),selection(:,2),20,'k','filled'); hold on
polarscatter(pi,max(selection(:,2)),200,'r','filled','MarkerFaceAlpha',.6);
polarscatter(pi,max(selection(:,2)),200,'r');
set(gca,'FontSize', 16)

%% Plot Boundary feature vector

files = dir('*_metadata.mat');      
load(files(3).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
idx = find(noFrames >= 8);
        
j = 2;
idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)} + pi/2;            
selection = cellCoordinates{idx(j),idx2(1)};
selection(:,1) = wrapTo2Pi(selection(:,1)); 
[x,y] = pol2cart(selection(:,1),selection(:,2));
mask = [round(x) round(y)];
ind = sub2ind([200,200], (mask(:,2)*-2)+100, (mask(:,1)*2)+100);
canvas = false(200,200);
canvas(ind) = true;
canvas = imfill(canvas,'holes');

I = double(cat(3, canvas, canvas, canvas));
I(I==0) = 0.99;
I(I==1) = 0.2;
s = 0.000001;
sumInd = boundary(mask,s);
mask2 = mask(sumInd(1:end-1),:);
mask2 = (mask2*diag([2 -2]))+100;
I2 = insertMarker(I,mask2,'Size',5, 'Color','blue');

vecSize = 15; % floor((size(mask2,1)/10)-1)*10;
set = mask2;
disp(size(mask2,1))
disp(vecSize)
while size(set,1) ~= vecSize
    % measure euclidean distance between 2 adjacent points and add a new 
    % point between the points with max distance    
    if size(set,1) < vecSize
        v1 = set;
        v2 = [set(2:end,:);set(1,:)];
        distance = sqrt(sum((v1-v2).^2,2));
        [~,idx3] = max(distance);
        if idx3 == size(set,1)
            newPoint = (set(end,:)+set(1,:))./2;
            set = [set;newPoint];
        else
            newPoint = (set(idx3,:)+set(idx3+1,:))./2;
            set = [set(1:idx3,:);newPoint;set(idx3+1:end,:)];
        end
    % measure euclidean disance between 3 adjacent points sequence and delete 
    % the point in the middle of the shortest distance   
    elseif size(set,1) > vecSize
        v1 = set;
        v2 = [set(3:end,:);set(1:2,:)];
        distance = sqrt(sum((v1-v2).^2,2));
        [~,idx3] = min(distance);
        if idx3 == size(set,1)
            set(1,:) = [];
        else
            set(idx3+1,:) = [];
        end
    end 
end

mask3 = set;
I3 = insertMarker(I,mask3,'square','Size',5,'Color','blue');
figure;
imshowpair(I2,I3,'montage');

%% Plot rotation for Rho vector

files = dir('*_metadata.mat');      
load(files(2).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     
seqLength = 8;            
vecSize   = 16;            

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
idx = find(noFrames >= 8);
        
j = 19;
idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)} - pi/2;            
selection = contourCoordinates{idx(j),idx2(1)};
selection(:,1) = wrapTo2Pi(selection(:,1) + rotation); 
contour = sortrows(selection,1);
contour = [contour ; contour(1,:)]; 

vecSize   = 12;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
rhoDescriptor = [];
for n = 1:vecSize
    [~,idx3] = min(abs(selection(:,1)-refAngles(n)));
    rhoDescriptor(n) = selection(idx3,2);
end

rhoLine = [refAngles' rhoDescriptor'];
rhoLine = [rhoLine ; rhoLine(1,:)];

figure;
polarplot(contour(:,1),contour(:,2),'Color', [.8 .8 .8],'lineWidth',12); hold on
polarscatter(0,max(selection(:,2)),200,'r','filled','MarkerFaceAlpha',.6);
polarscatter(0,max(selection(:,2)),200,'r');
polarplot(rhoLine(:,1),rhoLine(:,2),'r','lineWidth',1); hold on
pax = gca;
pax.ThetaAxisUnits = 'radians';
rticks([0 10 30 40])
set(gca,'FontSize', 16)

%% Plot difficult Rho shapes

files = dir('*_metadata.mat');      
load(files(2).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     

vecSize   = 16;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
idx = find(noFrames >= 8);
        
JJ = [19 5 7 23 45 34 65 27 43 36 33];

for k = 1:numel(JJ)
    j = JJ(k);
    idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
    rotation = rotationUp{idx(j),idx2(1)} - pi/2;            
    selection = contourCoordinates{idx(j),idx2(1)};
    selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);

    rhoDescriptor = zeros(vecSize,1);
    for n = 1:vecSize
        [~,idx3] = min(abs(selection(:,1)-refAngles(n)));
        rhoDescriptor(n,1) = selection(idx3,2);
    end    
    figure;
    ax = polaraxes;       
    polarscatter(ax,selection(:,1), selection(:,2),'k','filled'); hold on;
    polarscatter(ax,refAngles,rhoDescriptor, 'SizeData',100);
    ax.FontSize = 18;
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    rticks([])

    for i = 1:size(rhoDescriptor,1)
        polarplot([0;refAngles(i)],[0;rhoDescriptor(i)],'Color','red');
    end
    hold off;
    
    rhoLine = [refAngles' rhoDescriptor];
    rhoLine = [rhoLine ; rhoLine(1,:)];
    figure;
    ax = polaraxes;       
    polarscatter(ax,selection(:,1), selection(:,2),'k','filled'); hold on
    polarplot(ax,rhoLine(:,1),rhoLine(:,2),'r','lineWidth',3); hold off
    ax.FontSize = 18;
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    rticks([])
end


%% Plot rho descriptors of increasing length

files = dir('*_metadata.mat');      
load(files(2).name,'metadata');                                  
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     
seqLength = 8;            
vecSize   = 16;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);

% detect sequence of length equal or greater than seqLength
noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
idx = find(noFrames >= 8);
        
j = 19;
idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
rotation = rotationUp{idx(j),idx2(1)} - pi/2;            
selection = contourCoordinates{idx(j),idx2(1)};
selection(:,1) = wrapTo2Pi(selection(:,1) + rotation); 
selection = sortrows(selection,1);
plotSelec = [selection ; selection(1,:)];

test   = [4 8 16 32 64];
for k = 1:numel(test)
    vecSize = test(k);
    refAngles = linspace(0,2*pi,vecSize+1);
    refAngles = refAngles(1:end-1);
    rhoDescriptor = zeros(vecSize,1);
    for n = 1:vecSize
        [~,idx3] = min(abs(selection(:,1)-refAngles(n)));
        rhoDescriptor(n,1) = selection(idx3,2);
    end    
    figure;
    ax = polaraxes;
    polarplot(ax,plotSelec(:,1),plotSelec(:,2), 'LineWidth',3); hold on;
    polarscatter(ax,refAngles,rhoDescriptor, 'SizeData',150);
    ax.FontSize = 18;
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    rticks([])
    for i = 1:size(rhoDescriptor,1)
        polarplot(ax,[0;refAngles(i)],[0;rhoDescriptor(i)],'Color','red');
    end
    hold off;
        if vecSize==16
        rhoLine = [refAngles' rhoDescriptor];
        rhoLine = [rhoLine ; rhoLine(1,:)];
        figure;
        ax = polaraxes;
        polarscatter(ax,selection(:,1), selection(:,2),'k','filled'); hold on
        polarplot(ax,rhoLine(:,1),rhoLine(:,2),'r','lineWidth',1); hold off
        ax.FontSize = 18;
        pax = gca;
        pax.ThetaAxisUnits = 'radians';
        rticks([])
        end
end

%% Histograms for data augmentation

load('singleCellArea.mat','singleCellArea');
cellAreas = sort(cell2mat(singleCellArea'));
figure; 
histogram(cellAreas,'FaceColor',[.7 .7 .7], 'EdgeColor', [.85 .85 .85]); hold on
limits = quantile(cellAreas,[0.01 0.99]);  
xlim([limits(1) limits(2)])
xlabel('Single Cell Area (pixels)')
ylabel('Count')
set(gca,'fontsize',16)
hold on;
lines = quantile(cellAreas,[0.4 0.5 0.6 0.7]);  

yl = ylim;
for idx = 1 : numel(lines)
    plot([lines(idx) lines(idx)], [0 yl(2)], 'LineWidth',2);
end
hold off;

%% data augmentation

load('TrainingAreas.mat','TrainingAreas');   
load('RhoDescriptors.mat','RhoDescriptors');   
load('RhoDescriptors04.mat','RhoDescriptors04');   
load('RhoDescriptors05.mat','RhoDescriptors05');   
load('RhoDescriptors06.mat','RhoDescriptors06');   
load('RhoDescriptors07.mat','RhoDescriptors07');   
load('singleCellArea.mat','singleCellArea');

vecSize   = 16;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);

cellAreas = sort(cell2mat(singleCellArea'));
QQ = quantile(cellAreas,[0.4 0.5 0.6 0.7]);  

sample = [20 40 60 80 100];
for k = 1:numel(sample)
        i = sample(k);
    rhoLine = [refAngles' RhoDescriptors{i,1}];
    rhoLine = [rhoLine ; rhoLine(1,:)];
    rhoLine04 = [refAngles' RhoDescriptors04{i,1}*0.85];
    rhoLine04 = [rhoLine04 ; rhoLine04(1,:)];
    rhoLine05 = [refAngles' RhoDescriptors05{i,1}*0.95];
    rhoLine05 = [rhoLine05 ; rhoLine05(1,:)];
    rhoLine06 = [refAngles' RhoDescriptors06{i,1}*1.05];
    rhoLine06 = [rhoLine06 ; rhoLine06(1,:)];
    rhoLine07 = [refAngles' RhoDescriptors07{i,1}*1.15];
    rhoLine07 = [rhoLine07 ; rhoLine07(1,:)];

    figure;
    polarplot(rhoLine(:,1),rhoLine(:,2),':k','lineWidth',2); hold on
    polarplot(rhoLine04(:,1),rhoLine04(:,2),'lineWidth',2); hold on
    polarplot(rhoLine05(:,1),rhoLine05(:,2),'lineWidth',2); 
    polarplot(rhoLine06(:,1),rhoLine06(:,2),'lineWidth',2); 
    polarplot(rhoLine07(:,1),rhoLine07(:,2),'lineWidth',2);
    rlim([0 max(RhoDescriptors07{i,1}*1.2)])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    set(gca,'fontsize',16)
    rticks([])
    legend('original','0.4 Quantile', '0.5 Quantile', '0.6 Quantile', '0.7 Quantile');
end

%% Plot dynamic images

load('dynamicImages.mat','dynamicImages');     
% load('dynamicImages04.mat','dynamicImages04');     
% load('dynamicImages05.mat','dynamicImages05');     
% load('dynamicImages06.mat','dynamicImages06');     
% load('dynamicImages07.mat','dynamicImages07');     

No = [10 50 2350 2370];
for i = 1:numel(No)
    n = No(i);
    canvas = cat(3, dynamicImages{n,1}*0.8, dynamicImages{n,1}*0.8, dynamicImages{n,1}*0.8);
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




