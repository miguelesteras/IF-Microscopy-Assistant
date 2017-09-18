%% Phase 4.1. Transform cell mask into Fourier descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell mask into Fourier descriptors.
%   Fourier descriptors as first described by Kuhl and Giardina in 
%   "Elliptic Fourier features of a closed contour" 
%   Computer Graphics and Image Processing 18:236-258 1982
%   
%   This code is a modification of: 
%   http://geekstack.net/index.php?page=fourier_descriptor
%   Copyright © Tobias Pohlen 2017
%   A basic implementation of the descriptor presented in:
%   O. D. Trier, A. K. Jain, T. Taxt, Feature Extraction Methods for 
%   Character Recognition - A Survey, Pattern Recognition, Vol. 29, No. 4, 
%   pp. 641 - 662 (1996).
%   ======================================================================

close all; clear;
% seqLength = Number of frames in learning sequence (x input + 1 target).
% C determines the number of coefficients (and the description vector
% size). P determines the number of sampling points used for the fourier
% transform (in percentage of total number of points). The transformation
% is set not to be rotation invariant.
seqLength = 8;              
CC = [4 8 10];                 
P = 0.5;          
rotationInva = false;

for p = 1:numel(CC)
    C = CC(p);

    fourierDescriptor = [];
    files = dir('*_metadata.mat');      
    num_files = length(files);
    for i = 1:num_files
        load(files(i).name,'metadata');                                  
        load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');
        load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 

        % detect sequence of length equal or greater than seqLength
        noFrames = sum(double(~cellfun(@isempty,cellCoordinates)),2);
        idx = find(noFrames >= seqLength);
        fouriertemp = cell(100,seqLength+1); 
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
                    x = round(x)+metadata.maxRadious+1; 
                    y = round(y*-1)+metadata.maxRadious+1;
                    cellIdx = sub2ind(size(canvas), y, x);
                    canvas(cellIdx) = true;
                    canvas = imfill(canvas,'holes');
                    [a,b,c,d,~] = ellipticFourierDescriptor(canvas, C, P, rotationInva);
                    fouriertemp{count,m} = [a;b;c;d];
                end
                fouriertemp{count,end} = canvas;
                count = count +1;
            end
        end
        fourierDescriptor = [fourierDescriptor;fouriertemp];
    end
    save(strcat('fourierDescriptorC',num2str(C),'.mat'),'fourierDescriptor');
end
