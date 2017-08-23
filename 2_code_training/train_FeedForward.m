%% Phase 5.4 Train an FeedForward Network from Fourier Descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a FeedForward network using the fourier descriptors 
%   representation of cell shapes. 
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
X = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),inputData,'UniformOutput',false)))';
Y = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),targetData,'UniformOutput',false)))';

%% 4. Training feedforward network

hiddenSizes = [500 200];                          % neurons in hidden layer
ffNet = feedforwardnet(hiddenSizes);        % create feedforward network
view(ffNet)

% Fine-tunning
ffNet.trainFcn                    = 'trainoss';     % One-step secant backpropagation
ffNet.divideFcn                   = 'dividerand';   % data division random
ffNet.divideParam .trainRatio     = 0.9;              % all data used for training
ffNet.divideParam.valRatio        = 0.1;
ffNet.divideParam.testRatio       = 0;
ffNet.trainParam.epochs           = 1000;           % maximum number of epochs to train (for early stopping)
ffNet.trainParam.showCommandLine  = false;          % do not generate command-line output
ffNet.trainParam.showWindow       = true;           % show training GUI
ffNet.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)
ffNet.trainParam.show             = 1;              % Epochs between displays
tic;                                                % start timer

[ffNet,FFrecord] = train(ffNet,X,Y);
sprintf('Training time: %0.0f %s',toc,'s')          % duration of training

