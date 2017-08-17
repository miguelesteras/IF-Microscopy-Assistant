%% Phase 4.1.1. Transform cell contours into fourier descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. transform the binary image of a cell contour into a vector 
%   containing the fourier descriptors.
%   ======================================================================

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_dynamicImages.mat'),'dynamicImages');
    load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');

   
    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordenates)),2);
    idx = find(noFrames >= seqLength);
    BoundaryDataSet = cell(100,seqLength); 
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordenates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            % the rotation ([rotation to north] + pi/2) makes further 
            % pixel from center of image face west. 
            % This way the boundary descriptor starts from that point and
            % goes around the shape clockwise.
            rotation = rotationUp{idx(j),idx2(k)} + pi/2;
            for m = 1:seqLength                
                selection = cellCoordenates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                mask = [round(x) round(y)];
                
                % Start of Binary Search algorithm implementation for
                % values of the shrinking factor 's'. target size vecSize
                % + 1 because start and end value in boundary vector are
                % the same (and will be removed).
                from = 0; to = 1; 
                for n = 1:5
                    s = (from + to)/2;
                    sumInd = boundary(mask,s);
                    if numel(sumInd) == vecSize+1
                        break
                    elseif from == to
                        break
                    elseif numel(sumInd) < vecSize+1
                        from = s;
                    elseif numel(sumInd) > vecSize+1
                        to = s;
                    end
                end
                % End of Binary Search algorithm implementation
                sumInd = boundary(mask,s);

                hull = mask(sumInd(1:end-1),:);    
                BoundaryDataSet{count,m} = hull;
            end           
            count = count+1;
        end               
    end
    %% standarize feature vector size
    
    BoundaryDataSetSt = cell(size(BoundaryDataSet));
    for k = 1:numel(BoundaryDataSet)
        set = BoundaryDataSet{k};
        while size(set,1) ~= vecSize
            % measure euclidean disance between 2 adjacent points and add a new 
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
        BoundaryDataSetSt{k} = set;
    end
    
    save(strcat(metadata.name,'_BoundaryDataSet.mat'),'BoundaryDataSet');
    save(strcat(metadata.name,'_BoundaryDataSetSt.mat'),'BoundaryDataSetSt');     
    
%     clearvars -except files num_files i
end
