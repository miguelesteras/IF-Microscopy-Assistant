%% Phase 5.7. Train KNN with Generic Fourier descriptors (GFD)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell contour into GFD. Code by Frederik Kratzert
%   https://uk.mathworks.com/matlabcentral/fileexchange/52643-fd-=-gfd-bw-m-n--implementation-of-the-generic-fourier-descriptors
%   Binary image centered using centerobject funtion by Frederik Kratzert. 
%   ======================================================================

% seqLength = Number of frames in learning sequence (x input + 1 target).
% M determines the number of coefficients (and the description vector
% size). N determines the number of sampling points used for the fourier
% transform (in percentage of total number of points). The transformation
% is set not to be rotation invariant.
seqLength = 8;              
M = 9;          % radial frequency                 
N = 10;          % angular frequency    
GFDDescriptors = [];

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
    idx = find(noFrames >= seqLength);
    GFDtemp = cell(100,seqLength+1); 
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
                colSub = round(x)+metadata.maxRadious+1; 
                rowSub = round(y)+metadata.maxRadious+1;
                cellIdx = sub2ind(size(canvas), rowSub, colSub);
                canvas(cellIdx) = true;
                canvas = bwareaopen(bwmorph(imfill(canvas,'holes'),'bridge'),10);
                canvas = centerobject(canvas);
                descriptors = gfd(canvas,M,N);   
                GFDtemp{count,m} = descriptors(1:end-1);                
            end
            % the last column stores the binary image of the target frame
            GFDtemp{count,end} = canvas;
            count = count +1;
        end
    end
    GFDDescriptors = [GFDDescriptors;GFDtemp];
    clearvars -except files num_files i GFDDescriptor seqLength M N
end

save('GFDDescriptors.mat','GFDDescriptors');     

%% Create a dictionary of GFD descriptors using all single cells





