%%  Create Cell Cluster Sequences
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. 
%   ======================================================================

summaryLength = [];
files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    
    % create cell array with the location of only single cell in full
    % view as defined to be between min and max cell area.
    fullCellLocation = cell(size(singleCellMask,1),1);   
    
    for j = 1:size(singleCellMask,1)
        selection = bwareafilt(singleCellMask{j},[metadata.minCellArea metadata.maxCellArea]); % select only cells in full view  
        pixels = struct2cell(regionprops(selection,'PixelIdxList'));
        fullCellLocation{j} = pixels;     
    end
    
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
                if isempty(simil) || value < size(fullCellLocation{k}{m},1)*0.8
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
    %figure; histogram(summaryLength)

    % total number of unique sequences (>2 frames)
    numSeq = size(lengths(lengths > 1),1);
    % total number of 8 frames sequences (potential training data size)
    trainSize = sum(lengths(lengths > 7)-7);  
    [metadata(:).trainingSize] = trainSize;
    [metadata(:).uniqueSequences] = numSeq;

    save(strcat(metadata.name,'_fullCellLocation.mat'),'fullCellLocation');
    save(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    save(strcat(metadata.name,'_metadata.mat'),'metadata');
end
