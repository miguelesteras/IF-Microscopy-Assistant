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
%   the polar coordinate system).
%   ======================================================================

% Number of frames in learning sequence (seqLength) and descriptor vector
% size (vecSize).
seqLength = 8;            
vecSize   = 16;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
RhoDescriptors = [];
RhoDescriptors04 = [];
RhoDescriptors05 = [];
RhoDescriptors06 = [];
RhoDescriptors07 = [];
TrainingAreas = [];
cellBodyTraining = [];
contourTraining = [];
% Calculate area values for cummulative quantiles, later used for training
% data synthesis
load('singleCellArea.mat','singleCellArea');
cellAreas = sort(cell2mat(singleCellArea'));
QQ = quantile(cellAreas,[0.4 0.5 0.6 0.7]);

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
    idx = find(noFrames >= seqLength);
    TrainingAreasTemp = zeros(100,seqLength);
    tempRhoDesc = cell(100,seqLength);
    tempCellCoor = cell(100,seqLength);
    tempContCoor = cell(100,seqLength);
    tempRho04 = cell(100,seqLength);    
    tempRho05 = cell(100,seqLength);
    tempRho06 = cell(100,seqLength);
    tempRho07 = cell(100,seqLength);

    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            % the rotation ([rotation to north] - pi/2) makes further pixel from center of image face East. 
            % The rho descriptors starts from that point and goes around the shape anticlockwise.
            rotation = rotationUp{idx(j),idx2(k)} - pi/2;
            for m = 1:seqLength     
                TrainingAreasTemp(count,m) = size(cellCoordinates{idx(j),idx2(k+m-1)},1);
                cellBody = cellCoordinates{idx(j),idx2(k+m-1)};
                cellBody(:,1) = wrapTo2Pi(cellBody(:,1) + rotation);
                selection = contourCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
                rhoDescriptor = zeros(vecSize,1);
                for n = 1:vecSize
                    [~,idx3] = min(abs(selection(:,1)-refAngles(n)));
                    rhoDescriptor(n,1) = selection(idx3,2);
                end                    
                tempRhoDesc{count,m} = rhoDescriptor;
                tempCellCoor{count,m} = cellBody;
                tempContCoor{count,m} = selection;
            end 
            % create synthetic Rho feature vectors, introducing small scale
            % variations
            selectionArea = size(cellCoordinates{idx(j),idx2(k)},1);
            FF = sqrt(QQ/selectionArea);
            tempRho04(count,:) = cellfun(@(x) x*FF(1), tempRhoDesc(count,:),'un',0);
            tempRho05(count,:) = cellfun(@(x) x*FF(2), tempRhoDesc(count,:),'un',0);
            tempRho06(count,:) = cellfun(@(x) x*FF(3), tempRhoDesc(count,:),'un',0);
            tempRho07(count,:) = cellfun(@(x) x*FF(4), tempRhoDesc(count,:),'un',0);
            count = count+1;
        end               
    end
    TrainingAreas = [TrainingAreas ; TrainingAreasTemp];
    RhoDescriptors = [RhoDescriptors ; tempRhoDesc];
    RhoDescriptors04 = [RhoDescriptors04 ; tempRho04];
    RhoDescriptors05 = [RhoDescriptors05 ; tempRho05];
    RhoDescriptors06 = [RhoDescriptors06 ; tempRho06];
    RhoDescriptors07 = [RhoDescriptors07 ; tempRho07];
    cellBodyTraining = [cellBodyTraining ; tempCellCoor];
    contourTraining = [contourTraining ; tempContCoor];
end
save('TrainingAreas.mat','TrainingAreas');   
save('RhoDescriptors.mat','RhoDescriptors');   
save('RhoDescriptors04.mat','RhoDescriptors04');   
save('RhoDescriptors05.mat','RhoDescriptors05');   
save('RhoDescriptors06.mat','RhoDescriptors06');   
save('RhoDescriptors07.mat','RhoDescriptors07');   
save('cellBodyTraining.mat','cellBodyTraining');
save('contourTraining.mat','contourTraining');
