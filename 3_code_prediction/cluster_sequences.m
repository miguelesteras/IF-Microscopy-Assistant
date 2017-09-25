%%  Create Cell Cluster Sequences
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. 
%   ======================================================================
clear; close all; clc

files = dir('*_metadata.mat');      
num_files = length(files);
clusterAreas = [];
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');
    tempCusterArea = cell(numel(clusterMask),1);
    
    for m = 1:numel(clusterMask)
        CC = bwconncomp(clusterMask{m});
        tempCusterArea{m} = cellfun('length',CC.PixelIdxList);   
    end 
    clusterAreas = [clusterAreas ; tempCusterArea];
end
save('clusterAreas.mat','clusterAreas');


summaryLength = [];
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');

    clusterLocation = cell(size(clusterMask,1),1);       
    for j = 1:size(clusterMask,1)
        selection = bwareafilt(clusterMask{j},[metadata.minCellArea metadata.maxCellArea*10]);  
        pixels = struct2cell(regionprops(selection,'PixelIdxList'));
        clusterLocation{j} = pixels;     
    end
    
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
                if isempty(simil) || value < size(clusterLocation{k}{m},1)*0.8
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
    summaryLength = [summaryLength ; lengths];
    %figure; histogram(summaryLength)

    save(strcat(metadata.name,'_clusterLocation.mat'),'clusterLocation');
    save(strcat(metadata.name,'_clusterSequences.mat'),'clusterSequences');   
    save(strcat(metadata.name,'_metadata.mat'),'metadata');
end
