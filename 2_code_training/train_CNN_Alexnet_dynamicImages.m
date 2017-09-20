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
inputData = dynamicImages;
% transforms grayscale images into rgb, to match alexnet structure
inputData = cellfun(@(x) (cat(3, x, x, x)),inputData,'UniformOutput',false);

load('RhoDescriptors.mat','RhoDescriptors');
load('RhoDescriptors04.mat','RhoDescriptors04');
load('RhoDescriptors05.mat','RhoDescriptors05');
load('RhoDescriptors06.mat','RhoDescriptors06');
load('RhoDescriptors07.mat','RhoDescriptors07');
outputData = [RhoDescriptors ; RhoDescriptors04 ; RhoDescriptors05 ;
            RhoDescriptors06 ; RhoDescriptors07];
    
inputSize = size(inputData{1,1});
outputsize = size(size(outputData{1,1},1));
samples = size(inputData,1);

% set apart a random set of samples for validation
valRatio = .2;                        
vIdx = randi([1 samples],1,round(samples*valRatio)); 
inputTrain = inputData;    
inputTrain(vIdx,:) = [];
randIdx = randperm(size(inputTrain,1));
Xtrain = inputTrain(randIdx,1);

outputTrain = outputData;    
outputTrain(vIdx,:) = [];
Ytrain = outputTrain(randIdx,end);

Xval = inputData(vIdx,1);
Yval = outputData(vIdx,end);

% Load a pretrained AlexNet network.
net = alexnet;
net.Layers

% Remove last 3 output layers, to be adapted to new data
layersTransfer = net.Layers(1:end-3);
net.layersTransfer

% Add new layers
layers = [...
    layersTransfer
    fullyConnectedLayer(100,'WeightLearnRateFactor',20,...
    'BiasLearnRateFactor',20,'Name','fc8')
    fullyConnectedLayer(outputsize,'WeightLearnRateFactor',20,...
    'BiasLearnRateFactor',20,'Name','fc9')
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

% predict
YtrainPred = predict(netTrain,Xtrain);
YvalPred = predict(netTrain,Xval);
YpreFrame = outputData(vIdx,end-1);
netPerfor(m) = sqrt(mse(gsubtract(Yval,YvalPred)));
basePerfor(m) = sqrt(mse(gsubtract(Yval,YpreFrame)));
    
Alex_performance = mat2str([netPerfor basePerfor]);
save('perf_Alex.mat','Alex_performance')
save('net_Alex.mat','trainNetwork')

