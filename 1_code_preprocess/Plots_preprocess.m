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
Areas = [];
files = dir('*_metadata.mat');
num_files = length(files);
for k = 1:num_files
    load(files(k).name,'metadata');    
    load(strcat(metadata.name,'_singleCellArea.mat'),'singleCellArea'); 
    Areas = [Areas ; singleCellArea];
end

count = floor(numel(Areas)/10);
cellAreas = sort(cell2mat(Areas'));         
figure; histogram(cellAreas,count); 
hold on;
lines = quantile(cellAreas,[0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9]);  
plot([median(cellAreas) median(cellAreas)], [0 330], 'k', 'LineWidth',2);

for idx = 1 : numel(lines)
    plot([lines(idx) lines(idx)], [0 330], 'r', 'LineWidth',1);
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
ylabel('Total No. of Sequences (>=2 frame)')
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

figure; histogram(summaryLength,90); hold on
yl = ylim;
plot([8 8], [0 yl(2)], 'r', 'LineWidth',2);
ylabel('Count')
xlabel('Unique Sequence Length (frames)')
set(gca,'fontsize',16)
