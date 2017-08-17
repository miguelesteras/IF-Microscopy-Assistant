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
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
                
files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_contourCoordenates.mat'),'contourCoordenates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,contourCoordenates)),2);
    idx = find(noFrames >= seqLength);
    RhoDescriptors = cell(100,seqLength); 
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
                rhoDescriptor = zeros(vecSize,1);
                for n = 1:vecSize
                    [~,idx3] = min(abs(selection-refAngles(n)));
                    rhoDescriptor(n) = selection(idx3,2);
                end                    
                RhoDescriptors{count,m} = rhoDescriptor;
            end           
            count = count+1;
        end               
    end    
    save(strcat(metadata.name,'_RhoDescriptors.mat'),'RhoDescriptors');
    
    %clearvars -except files num_files i
end

%%
% 
% figure
% selection = contourCoordenates{3,1};
% rotation = rotationUp{3,1};
% 
% polarplot(selection(:,1),selection(:,2)); hold on
% % selection2 = [wrapTo2Pi(selection(:,1) + rotation) selection(:,2)];
% polarplot(selection2(:,1),selection2(:,2));
% selection3 = [wrapTo2Pi(selection(:,1) + rotation - pi/2)  selection(:,2)];
% polarplot(selection3(:,1),selection3(:,2));
