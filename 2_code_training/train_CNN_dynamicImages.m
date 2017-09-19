%% Phase 5.7 Train CNN for dynamic Images (based on AlexNet in Matlab)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Use pre-trained CNN to extract featuers from dynamic images and 
%   predict with them rho descriptors .
%   ======================================================================
clear; close all; clc;
rng('default')

load('dynamicImages_Alexnet.mat','dynamicImages');
load('RhoDescriptors.mat','RhoDescriptors');
inputData = dynamicImages;

% transforms grayscale images into rgb, to match alexnet structure
inputData = cellfun(@(x) (cat(3, x, x, x)),inputData,'UniformOutput',false);

outputData = RhoDescriptors;
inputSize = size(dynamicImages{1,1});
outputsize = numel(RhoDescriptors{1,1});
samples = size(data,1);
[tIdx, vIdx] = dividerand(samples, .85, .15);

Xtrain = inputData(tIdx,1);
Ytrain = outputData(tIdx,end);
Xval = inputData(vIdx,1);
Yval = outputData(vIdx,end);

% Load a pretrained AlexNet network.
net = alexnet;
net.Layers

% Remove last 3 output layers, to be adapted to new data
layersTransfer = net.Layers(1:end-3);
net.Layers

% Add new layers
layers = [...
    layersTransfer
    fullyConnectedLayer(100,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    fullyConnectedLayer(outputsize,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    ];

% training options
options = trainingOptions('sgdm',...        % stochastic gradient descent with momentum
    'MiniBatchSize',5,...                   % batch size (no. pictures)
    'MaxEpochs',10,...                      % max. epochs for early stopping
    'InitialLearnRate',0.001,...            % initial learning rate
    'LearnRateSchedule', 'piecewise', ...   % the learning rate decreases
    'LearnRateDropFactor', 0.1, ...         % by factor of 0.1
    'LearnRateDropPeriod', 1);              % every 8 epochs

% re-train Alexnet
netTrain = trainNetwork(Xtrain, Ytrain, layers, options);

layer = 'fc9';      % last fully connected layer form AlexNet (1000 neurons)
trainFeatures = activations(netTrain,Xtrain,layer);     % extract features for training data
testFeatures = activations(netTrain,Xtest,layer);       % ectract features for test data



