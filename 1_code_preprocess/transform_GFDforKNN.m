%% Phase 4.7. Transform cell mask coordinates into Generic Fourier descriptors (GFD)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell coordinates into GFD. Code by Frederik Kratzert
%   https://uk.mathworks.com/matlabcentral/fileexchange/52643-fd-=-gfd-bw-m-n--implementation-of-the-generic-fourier-descriptors
%   Binary image centered using centerobject funtion by Frederik Kratzert. 
%   ======================================================================

% seqLength = Number of frames in learning sequence (x input + 1 target).
% M is the radial frequency. N is the angular frequency.
seqLength = 8;              
M = 9;                           
N = 10;   
GFDDescriptors = [];
GFDDescriptorsScaled = [];
scaleFactor = [linspace(1,2,seqLength-1) 1];

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    
    % GFdtemp stores descriptor vectors as calculated by the function GFD.
    % In GFDtemp2 the values in the vectors are scaled up in funtion of
    % time with a scale factor in the range [1 2]. 1 for the first input 
    % frame in the series and 2 for the last input frame. The output frame
    % is not scaled up. Hence a KNN classifier using a distance measure
    % like 'euclidean' will weight higher frames closer to the output.
    GFDtemp = cell(100,seqLength+1);
    GFDtemp2 = cell(100,seqLength+1); 

    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,cellCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};            
            for m = 1:seqLength
                canvas = false((metadata.maxRadious*2)+1, (metadata.maxRadious*2)+1);
                selection = cellCoordinates{idx(j),idx2(k+m-1)};
                selection(:,1) = wrapTo2Pi(selection(:,1) + rotation);          
                [x,y] = pol2cart(selection(:,1),selection(:,2));
                % including a change in the coordenates system, from origin
                % [0,0] being center of image to [0,0] being top left image
                % corner and postive y axis = rows (hence 'y*-1' is needed)
                colSub = round(x)+metadata.maxRadious+1; 
                rowSub = round(y*-1)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                canvas = bwareaopen(imfill(bwmorph(canvas,'bridge'),'holes'),10);                
                canvas = centerobject(canvas);
                descriptors = gfd(canvas,M,N);   
                GFDtemp{count,m} = descriptors(1:end-1);
                GFDtemp2{count,m} = descriptors(1:end-1)*scaleFactor(m);
            end
            % the last column stores the binary image of the target frame
            GFDtemp2{count,end} = canvas;
            GFDtemp{count,end} = canvas;
            count = count +1;
        end
    end
    GFDDescriptors = [GFDDescriptors;GFDtemp];
    GFDDescriptorsScaled = [GFDDescriptorsScaled;GFDtemp2];
end

save('GFDDescriptors.mat','GFDDescriptors');     
save('GFDDescriptorsScaled.mat','GFDDescriptorsScaled');     