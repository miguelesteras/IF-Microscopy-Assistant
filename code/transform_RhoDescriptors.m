%% Phase 4.4. Transform cell contours into a rho (distance from center) descriptor vector
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell contour into a descriptor vecctor containing the rho
%   value for a given number of contour pixels (distance from center, in 
%   the polar coordenate system).
%   ======================================================================

%% Build data set of single cells

seqLength = 8;             % Number of frames in learning sequence (x input + 1 target)
vecSize   = 32;            % Descriptor vector size

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_cellSequences.mat'),'cellSequences');   
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_centerMass.mat'),'centerMass');  
    load(strcat(metadata.name,'_cellCoordenates.mat'),'cellCoordenates');     
    load(strcat(metadata.name,'_contourCoordenates.mat'),'contourCoordenates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,contourCoordenates)),2);
    idx = find(noFrames >= seqLength);
    BoundaryDataSet = cell(100,seqLength); 
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,contourCoordenates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            % the rotation ([rotation to north] - pi/2) makes further pixel from center of image face East. 
            % The rho descriptors starts from that point and goes around the shape anticlockwise.
            rotation = rotationUp{idx(j),idx2(k)} - pi/2;
            for m = 1:seqLength                
                selection = contourCoordenates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
                
                % extract contour pixels only
                % boundary(x,y,s)
                %
                % 
                refAngles = linspace(0,2*pi,vecSize+1);
                refAngles = refAngles(1:end-1);
                rhoDescriptor = zeros(vecSize,1);
                for n = 1:numel(refAngles)
                    [~,idx3] = min(abs(selection-refAngle(n)));
                    rhoDescriptor(n) = selection(idx3,2);
                
                
                
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
    %%
    % standarize feature vector size
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
    
    %clearvars -except files num_files i
end

%%

figure
selection = contourCoordenates{3,1};
rotation = rotationUp{3,1};

polarplot(selection(:,1),selection(:,2)); hold on
selection2 = [wrapTo2Pi(selection(:,1) + rotation) selection(:,2)];
polarplot(selection2(:,1),selection2(:,2));
selection3 = [wrapTo2Pi(selection(:,1) + rotation - pi/2)  selection(:,2)];
polarplot(selection3(:,1),selection3(:,2));


%%

figure 
ind = sub2ind([200,200], y+100, x+100);
canvas = false(200,200);
canvas(ind) = true;