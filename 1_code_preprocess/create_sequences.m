%% Phase 2. Create Cell Sequences
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. create a data set with sequencial frames of single cells 
%   ======================================================================

%% Build data set of single cells

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_singleCellMask.mat'),'singleCellMask');
    
    fullCellMask = cell(size(singleCellMask,1),1);      % initialize mask for full view single cells 
    fullCellLocation = cell(size(singleCellMask,1),1);  % initialize record for pixel location of full view cells 
    
    for j = 1:size(singleCellMask,1)
        selection = bwareafilt(singleCellMask{j},[metadata.minCellArea metadata.maxCellArea]); % select only cells in full view  
        fullCellMask{j} = selection;
        pixels = struct2cell(regionprops(selection,'PixelIdxList'));
        fullCellLocation{j} = pixels;     
    end
    
    save(strcat(metadata.name,'_fullCellMask.mat'),'fullCellMask');     % save files on disk
    save(strcat(metadata.name,'_fullCellLocation.mat'),'fullCellLocation');
    
    %%
    cellSequences = cell(1,size(singleCellMask,1)); % cell array to store cell sequences
    seqCount = 1;
    
    for m = 1:size(fullCellLocation{1},2)   % detect cells and start sequence in first frame
        object = fullCellLocation{1}{m};
        cellSequences{seqCount,1} = object;
        seqCount = seqCount + 1;
    end
            
    for k = 2:size(fullCellLocation,1)  % from second frame onwards
        if isempty(fullCellLocation{k})
            continue
        else
            preSeq = find(~cellfun(@isempty,cellSequences(:,k-1)));
            prior = cellSequences(preSeq,k-1);

            for m = 1:size(fullCellLocation{k},2)   % for every cell in frame
                object = fullCellLocation{k}(m);
                object = repmat(object,[size(prior,1) 1]);
                simil = cellfun('length',cellfun(@intersect, prior, object, 'UniformOutput', false));
                [value, idx] = max(simil);          % most similar cell in previous frame
                
                % if cell does not match previous => 60% location 
                % similarity or previous frame is empty start new sequence:
                if isempty(simil) || value < size(fullCellLocation{k}{m},1)*0.6
                    cellSequences{seqCount,k} = fullCellLocation{k}{m};
                    seqCount = seqCount + 1;
                % else add cell to a previous sequence
                else    
                    cellSequences{preSeq(idx),k} = fullCellLocation{k}{m};                      
                end
            end
        end
    end
    binarymap = ~cellfun(@isempty,cellSequences);
    figure; imshow(binarymap);                      % show binary map of sequences
    seqLength = sum(uint8(binarymap),2);
    %figure; histogram(seqLength)                    % show frequency histogram of sequences
    numSeq = size(seqLength(seqLength > 4),1);      % number of unique sequences > 4 frames
    trainSize = sum(seqLength(seqLength > 4)-4);    % total number of 5 frames sequences

    [metadata(:).uniqueSequences] = numSeq;     % add info to metadata file
    [metadata(:).trainingSize] = trainSize;

    save(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   % save files on disk
    save(strcat(metadata.name,'_metadata.mat'),'metadata');

    clearvars -except files num_files i
end
