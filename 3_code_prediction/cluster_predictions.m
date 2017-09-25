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
for i = 1:num_files
    load(files(i).name,'metadata');
    load(strcat(metadata.name,'_clusterMask.mat'),'clusterMask');
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');
    load(strcat(metadata.name,'_clusterLocation.mat'),'clusterLocation');
    load(strcat(metadata.name,'_movie.mat'),'movie');
    % cellPredictions contains the cell sequences and the prediction of
    % single cell shape part of cell clusters.
    cellPredictions = cellSequences;
    
    for j = 1:size(cellPredictions,1)
        % next frame to last single cell in sequence 
        idx = find(~cellfun(@isempty,cellPredictions(j,:)), 1, 'last' )+1;
        
        % if there are no clusters in the next frame, step to next sequence
        if isempty(clusterLocation{idx})
            continue
        % else, find if single cell is part of a cluster
        else
            clusters = clusterLocation{idx};
            cell = cellPredictions(j,idx-1);
            cell = repmat(cell,[1 size(clusters,2)]);
            simil = cellfun('length',cellfun(@intersect, cell, clusters, 'UniformOutput', false));
            [value, idx2] = max(simil);
            % if single cell is not part of a clusters, step to next sequence
            if isempty(simil) || value < size(cellPredictions{j,idx-1},1)*0.8
                continue
            else
                while true
                    frame = movie{idx};
                    % predict single cell shape in next frame
                    cellPredictions{j,idx} = predictShape(cell,...
                                                          clusters(idx2),...
                                                          frame);                                                                                          
                    idx = idx+1;
                    % if there are no clusters in the next frame, stop while loop
                    if isempty(clusterLocation{idx})
                        break
                    else
                        clusters = clusterLocation{idx};
                        cell = cellPredictions(j,idx-1);
                        cell = repmat(cell,[1 size(clusters,2)]);
                        simil = cellfun('length',cellfun(@intersect, cell, clusters, 'UniformOutput', false));
                        [value, idx2] = max(simil);
                        % if single cell is not part of a clusters, stop wile loop
                        if isempty(simil) || value < size(cellPredictions{j,idx-1},1)*0.8
                            break 
                        end
                    end                               
                end
            end
        end
    end
end

        
        
               