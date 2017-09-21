%% Phase 5.7 Build and Train CNN for dynamic Images
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Create a new CNN and train it to extract featuers from dynamic images and 
%   predict with them rho descriptors .
%   ======================================================================
clear; close all; clc;
rng('default')

load('dynamicImages.mat','dynamicImages');
load('dynamicImages04.mat','dynamicImages04');
load('dynamicImages05.mat','dynamicImages05');
load('dynamicImages06.mat','dynamicImages06');
load('dynamicImages07.mat','dynamicImages07');
inputData = [dynamicImages ; dynamicImages04 ; dynamicImages05 ;
            dynamicImages06 ; dynamicImages07];

load('RhoDescriptors.mat','RhoDescriptors');
load('RhoDescriptors04.mat','RhoDescriptors04');
load('RhoDescriptors05.mat','RhoDescriptors05');
load('RhoDescriptors06.mat','RhoDescriptors06');
load('RhoDescriptors07.mat','RhoDescriptors07');
outputData = [RhoDescriptors ; RhoDescriptors04 ; RhoDescriptors05 ;
            RhoDescriptors06 ; RhoDescriptors07];

inputSize = size(inputData{1,1});
outputsize = size(outputData{1,1},1);
samples = size(inputData,1);

% set apart a random set of samples for validation
valRatio = .2;                        
vIdx = randi([1 samples],1,round(samples*valRatio)); 
inputTrain = inputData;    
inputTrain(vIdx,:) = [];
randIdx = randperm(size(inputTrain,1));
Xtrain = inputTrain(randIdx,1);
% convet intput from cell array to 4-D matrix
Xtrain = cat(4, Xtrain{:});

outputTrain = outputData;    
outputTrain(vIdx,:) = [];
Ytrain = outputTrain(randIdx,end);

Xval = inputData(vIdx,1);
% convet from cell array to 4-D matrix
Xval = cat(4, Xval{:});
Yval = outputData(vIdx,end);

% 1 convolution

% define layers
layers = [imageInputLayer([inputSize 1])
          convolution2dLayer(10,20)
          reluLayer
          maxPooling2dLayer(3,'Stride',3)
          fullyConnectedLayer(100,'Name','fc1')
          fullyConnectedLayer(16,'Name','fc2')
          regressionLayer];      
      
% training options
options = trainingOptions('sgdm',...        % stochastic gradient descent with momentum
    'MiniBatchSize',3,...                   % batch size (no. pictures)
    'MaxEpochs',100,...                      % max. epochs for early stopping
    'InitialLearnRate',0.001,...            % initial learning rate
    'LearnRateSchedule', 'piecewise', ...   % the learning rate decreases
    'LearnRateDropFactor', 0.1, ...         % by factor of 0.1
    'LearnRateDropPeriod', 1);              % every 8 epochs

% Train CNN
[net,record] = trainNetwork(Xtrain, Ytrain, layers, options);

% predict
YtrainPred = predict(net,Xtrain);
YvalPred = predict(net,Xval);
YpreFrame = outputData(vIdx,end-1);
netPerfor(m) = sqrt(mse(gsubtract(Yval,YvalPred)));
basePerfor(m) = sqrt(mse(gsubtract(Yval,YpreFrame)));
    
CNN_performance = mat2str([netPerfor basePerfor]);
save('perf_CNN1.mat','CNN_performance')
save('record_CNN1.mat','record')
save('net_CNN1.mat','net')


% 2 concolutions

% define layers
layers = [imageInputLayer([inputSize 1])
          convolution2dLayer(10,20)
          reluLayer
          maxPooling2dLayer(3,'Stride',3)
          convolution2dLayer(10,20)
          reluLayer
          maxPooling2dLayer(2,'Stride',2)
          fullyConnectedLayer(100,'Name','fc1')
          fullyConnectedLayer(16,'Name','fc2')
          regressionLayer];
          
% training options
options = trainingOptions('sgdm',...        % stochastic gradient descent with momentum
    'MiniBatchSize',3,...                   % batch size (no. pictures)
    'MaxEpochs',100,...                      % max. epochs for early stopping
    'InitialLearnRate',0.001,...            % initial learning rate
    'LearnRateSchedule', 'piecewise', ...   % the learning rate decreases
    'LearnRateDropFactor', 0.1, ...         % by factor of 0.1
    'LearnRateDropPeriod', 1);              % every 8 epochs

% Train CNN
[net,record] = trainNetwork(Xtrain, Ytrain, layers, options);

% predict
YtrainPred = predict(net,Xtrain);
YvalPred = predict(net,Xval);
YpreFrame = outputData(vIdx,end-1);
netPerfor(m) = sqrt(mse(gsubtract(Yval,YvalPred)));
basePerfor(m) = sqrt(mse(gsubtract(Yval,YpreFrame)));
    
CNN_performance = mat2str([netPerfor basePerfor]);
save('perf_CNN2.mat','CNN_performance')
save('record_CNN2.mat','record')
save('net_CNN2.mat','net')
