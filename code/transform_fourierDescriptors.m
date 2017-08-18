%% Phase 4.1. Transform cell silhouettes into Fourier descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Extract sequences of frames (f1, f2, ... , fx), from which
%   input sequence = (f1, f2, ... , fx-1), and target sequence = fx
%   2. Transform cell silhouette into Fourier descriptors.
%   Fourier descriptors as first described by Kuhl and Giardina in 
%   "Elliptic Fourier features of a closed contour" 
%   Computer Graphics and Image Processing 18:236-258 1982
%   This implementation has been created by David Thomas from the 
%   University of Melbourne, and can be found at:
%   https://uk.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors
%   Copyright (c) 2005, David Thomas
%   ======================================================================

%% Build data set of single cells

seqLength = 8;              % Number of frames in learning sequence (x input + 1 target)
NoHarmonics = 25;           % Number of Harmonics (interger greater than 0)
NormSize = true;            % normalize size of object
NormOrientation = false;    % normalize orientation of object

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

%% Fourier descriptor using code from:
% http://geekstack.net/index.php?page=fourier_descriptor
% Copyright � Tobias Pohlen 2017
% A basic implementation of the descriptor presented in:
% O. D. Trier, A. K. Jain, T. Taxt, Feature Extraction Methods for 
% Character Recognition - A Survey, Pattern Recognition, Vol. 29, No. 4, 
% pp. 641 - 662 (1996).

load(files(1).name,'metadata');                                  
load(strcat(metadata.name,'_cellCoordinates.mat'),'cellCoordinates'); 
load(strcat(metadata.name,'_rotationUp.mat'),'rotationUp');

input = cellCoordinates{3,1};
input(:,1) = input(:,1) + rotationUp{3,1};   
[x,y] = pol2cart(input(:,1),input(:,2));
selection = [round(x)+100 round(y)+100];         
ind = sub2ind([200,200], selection(:,2), selection(:,1));
canvas = false(200,200);
canvas(ind) = true;

% fourier descriptor
C = 20;                 % number of coefficients (description vector size)
P = 0.12;               % number of sampling points in percent
rotationInva = false;
[a,b,c,d,~] = ellipticFourierDescriptor(canvas, C, P, rotationInva);

% Image reconstruction form fourier descriptor. T = number of sampling 
% points for reconstruction.
xt = ellipticFourierFunction(a,b,T);
yt = ellipticFourierFunction(c,d,T);
plp = zeros(1,T);
for t = 1:T
    plp(t) = xt(t) + i*yt(t);
end
% plot reconstructed image with center [100,100] and rescale by scaleF.
scaleF = 50;
x = round(real(plp)*scaleF)';
y = round(imag(plp)*scaleF)';
ind2 = sub2ind([200,200], y+100, x+100);
canvas2 = false(200,200);
canvas2(ind2) = true;
figure
imshowpair(canvas,canvas2,'montage');
