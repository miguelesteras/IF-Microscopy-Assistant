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
%   ======================================================================

% Number of frames in learning sequence (x input + 1 target)
seqLength = 8; 

files = dir('*_metadata.mat');      
num_files = length(files);
for i = 1:num_files
    load(files(i).name,'metadata');                                  
    load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
    load(strcat(metadata.name,'_contourCoordinates.mat'),'contourCoordinates');     

    % detect sequence of length equal or greater than seqLength
    noFrames = sum(double(~cellfun(@isempty,contourCoordinates)),2);
    idx = find(noFrames >= seqLength);
    fourierDescriptor = cell(100,seqLength); 
    count = 1;
    
    for j = 1:size(idx,1)
        idx2 = find(~cellfun(@isempty,contourCoordinates(idx(j),:)));
        for k = 1:size(idx2,2) - (seqLength-1)
            rotation = rotationUp{idx(j),idx2(k)};
            
            for m = 1:seqLength                
                input = contourCoordinates{idx(j),idx2(k+m-1)};
                input(:,1) = input(:,1) + rotation;          % apply rotation
                [x,y] = pol2cart(input(:,1),input(:,2));
                contour = [round(x) round(y)];         
                descriptor = fEfourier(contour, NoHarmonics, NormSize, NormOrientation);
                fourierDescriptor{count,m} = descriptor;
            end        
        count = count +1;
        end
    end
    save(strcat(metadata.name,'_fourierDescriptor.mat'),'fourierDescriptor');     

    %clearvars -except files num_files i
end

%% grphs

ind = sub2ind([200,200], contour(:,2)+100, contour(:,1)+100);
canvas = false(200,200);
canvas(ind) = true;
% reverse descriptor
descriptor = fEfourier(contour, NoHarmonics, NormSize, NormOrientation);
image = rEfourier(descriptor, NoHarmonics, 100);
ind2 = sub2ind([200,200], (round(image(:,2)/10)+100), round((image(:,1)/10)+100));
canvas2 = false(200,200);
canvas2(ind2) = true;
figure
imshowpair(canvas,canvas2,'montage');


