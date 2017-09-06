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

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
    idx = find(noFrames >= seqLength);
    tempRhoDesc = cell(100,seqLength); 
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            % the rotation ([rotation to north] - pi/2) makes further pixel from center of image face East. 
            % The rho descriptors starts from that point and goes around the shape anticlockwise.
            rotation = rotationUp{idx(j),idx2(k)} - pi/2;
            for m = 1:seqLength                
                selection = contourCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
                rhoDescriptor = zeros(vecSize,1);
                for n = 1:vecSize
                    [~,idx3] = min(abs(selection(:,1)-refAngles(n)));
                    rhoDescriptor(n,1) = selection(idx3,2);
                end                    
                tempRhoDesc{count,m} = rhoDescriptor;
            end           
            count = count+1;
        end               
    end
    %save(strcat(metadata.name,'_RhoDescriptors.mat'),'tempRhoDesc');
    RhoDescriptors = [RhoDescriptors;tempRhoDesc];
end
save('RhoDescriptors.mat','RhoDescriptors');   

%% Plots
% rotation = -5.977;
% selection = contourCoordinates{idx(2),idx2(2)};
% selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);
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
    polarplot(plotSelec(:,1),plotSelec(:,2), 'LineWidth',3); hold on;
    polarscatter(refAngles,rhoDescriptor, 'SizeData',100);
    for i = 1:size(rhoDescriptor,1)
        polarplot([0;refAngles(i)],[0;rhoDescriptor(i)],'Color','red');
    end
    hold off;
end
