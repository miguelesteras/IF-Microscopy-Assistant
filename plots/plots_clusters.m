%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Plots for cluster Areas

load('clusterAreas.mat','clusterAreas');
clusterAreas = sort(cell2mat(singleCellArea'));
count = floor(numel(clusterAreas)/10);
figure; histogram(clusterAreas,count, 'FaceColor',[.7 .7 .7], 'EdgeColor', [.85 .85 .85]); 
limits = quantile(clusterAreas,[0.001 0.99]);  
xlim([limits(1) limits(2)])
xlabel('Cluster Area (pixels)')
ylabel('Count')
set(gca,'fontsize',16)
hold on;
lines = quantile(clusterAreas,[0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9]);  

yl = ylim;
plot([median(clusterAreas) median(clusterAreas)], [0 yl(2)], 'k', 'LineWidth',2);

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


%% Plots on cluster sequences

summaryLength = [];
files = dir('*_metadata.mat');      
load(files(2).name,'metadata');
load(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');
load(strcat(metadata.name,'_clusterLocation.mat'),'clusterLocation');

MeanLength = [];
NumSequence = [];
NumSeq8     = [];

X = [0.5 0.6 0.7 0.75 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 0.99];
for i = 1:numel(X)  
    S = X(i);

    % cell array to store cell sequences
    clusterSequences = cell(1,size(clusterMask,1)); 
    seqCount = 1;

    % detect cells and start sequence in first frame
    for m = 1:size(clusterLocation{1},2)   
        object = clusterLocation{1}{m};
        clusterSequences{seqCount,1} = object;
        seqCount = seqCount + 1;
    end

    % from second frame onwards, for every cell in frame, check if the cell
    % was present in previoius frame (based on >= similar location). If
    % true, add cell to previous sequence, if false, start a new sequence.
    for k = 2:size(clusterLocation,1)  
        if isempty(clusterLocation{k})
            continue
        else
            preSeq = find(~cellfun(@isempty,clusterSequences(:,k-1)));
            prior = clusterSequences(preSeq,k-1);
            for m = 1:size(clusterLocation{k},2)
                object = clusterLocation{k}(m);
                object = repmat(object,[size(prior,1) 1]);
                simil = cellfun('length',cellfun(@intersect, prior, object, 'UniformOutput', false));
                [value, idx] = max(simil);                          
                if isempty(simil) || value < size(clusterLocation{k}{m},1)*S
                    clusterSequences{seqCount,k} = clusterLocation{k}{m};
                    seqCount = seqCount + 1;
                else    
                    clusterSequences{preSeq(idx),k} = clusterLocation{k}{m};                      
                end
            end
        end
    end

    % Binary map of sequences, histogram and sequence lengths
    binarymap = ~cellfun(@isempty,clusterSequences);
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

% sequence length histogram
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

