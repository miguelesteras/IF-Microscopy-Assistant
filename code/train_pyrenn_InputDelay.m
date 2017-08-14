%% Phase 5.5 Train Recurrent Network with input delay (pyrenn)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a Recurrent network with input delay using the pyrenn
%   package and the fourier descriptors.
%   ======================================================================

rng('default')

% Create a joined dataSet from all video files
files = dir('*_metadata.mat');      
num_files = length(files);
inputData = [];   
targetData = [];
for i = 1:num_files
    load(files(i).name,'metadata');                               
    load(strcat(metadata.name,'_fourierInput.mat'),'fourierInput');     
    load(strcat(metadata.name,'_fourierTarget.mat'),'fourierTarget');     
    inputData = [inputData; fourierInput];
    targetData = [targetData; fourierTarget];
end

% transform input and output from cell array to column vectors
inputData = (cellfun(@(x) (reshape(x', [], 1)),inputData,'UniformOutput',false));
targetData = (cellfun(@(x) (reshape(x', [], 1)),targetData,'UniformOutput',false));
dataSet = [inputData targetData]; 

% Create network structure with delays in input only
inputSiz    = 100;
hidden1Siz  = 500;
hidden2Siz  = 500;
outputSiz   = 100;
dIn         = [0,1,2,3];  % Set of input delays of the neural network
dIntern     = [];  % Set of inernal delays of the neural network
dOut        = [];  % Set of output delays of the neural network

net = CreateNN([inputSiz hidden1Siz hidden2Siz outputSiz], dIn, dIntern, dOut);

%% Train network with training data, one input at the time.
repeat      = 2;    % training events with same input sequence
terError    = 1e-5; % target prediction error to terminate
epochs      = 500;  % number of epochs

for j = 1:epochs
    for k = 1:size(dataSet,1)
        X = cell2mat(dataSet(k,1:4));
        Y = cell2mat(dataSet(k,2:5));
        net = train_LM(X,Y,net,repeat,terError);
    end    
end

%% Validate
yHat = NNOut(X,net);




