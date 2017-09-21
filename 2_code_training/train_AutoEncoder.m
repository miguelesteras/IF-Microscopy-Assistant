%% Phase 5.3 Train an AutoEncoder from Summary Image
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a AutoEncoder network for summary image input, followed by a 
%   fully connected layer to form a output frame. 
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

% Compute a Gaussian pyramid reduction by one level
for j = 1:size(inputData,1)
    inputData{j,1} = impyramid(inputData{j,1}, 'reduce');    
end

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
valRatio = .15;                        
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

%% Train a autoencoder for dynamic images: First layer
enco_fnc = ['logsig' 'satlin'];
deco_fnc = ['logsig' 'satlin' 'purelin'];

hiddenSize1 = 1000;                         	
AE1 = trainAutoencoder(Xtrain,hiddenSize1, ...  % images in cell array format
    'MaxEpochs',500, ...                       % max no.epochs for early stopping
    'EncoderTransferFunction', 'logsig',...     % logistic sigmoid function
    'DecoderTransferFunction', 'logsig',...     % logistic sigmoid function
    'LossFunction', 'msesparse',...             % loss function used for training
    'TrainingAlgorithm', 'trainscg',...         % training algorithm, scaled conjugate gradient descent
    'L2WeightRegularization',0.004, ...         % coefficient for the L2 weight regularizer in the cost function
    'SparsityRegularization',4, ...             % controls the weighting of the sparsity regularizer
    'SparsityProportion',0.2, ...               % proportion of training examples which a neuron in the hidden layer should activate in response to.
    'ScaleData', false);                        % do not rescale input data

AE1features = encode(AE1,Xtrain);            % generate feautres from autoencoder1 hidden layer
save('net_AE1.mat','AE1')

%% Second layer

% use features extracted by autoencoder1 as input to autoencoder 2

hiddenSize2 = 500;
AE2 = trainAutoencoder(AE1features,hiddenSize2, ... 
    'MaxEpochs',500, ...
    'EncoderTransferFunction', 'logsig',...
    'DecoderTransferFunction', 'logsig',...
    'LossFunction', 'msesparse',...
    'TrainingAlgorithm', 'trainscg',...                 % training algorithm, scaled conjugate gradient descent
    'L2WeightRegularization',0.004, ...                 % coefficient for the L2 weight regularizer in the cost function
    'SparsityRegularization',4, ...
    'SparsityProportion',0.2, ...
    'ScaleData', false);                                % do not rescale input data

save('net_AE2.mat','AE2')

%% Train deep network with fully connected layer

% transform input and output from cell array to column vectors
X = cell2mat(cellfun(@(x) (reshape(x, [], 1)),Xtrain','UniformOutput',false));
Y = cell2mat(cellfun(@(x) (reshape(x, [], 1)),Ytrain','UniformOutput',false));

% train network stack
FC1 = fullyConnectedLayer(200);
FC2 = fullyConnectedLayer(100);
FC3 = fullyConnectedLayer(16);

network = stack(AE1,AE2,FC);
view(network)

network.inputs{1}.size          = 25281;
network.trainFcn                = 'trainoss';   % One-step secant backpropagation
network.trainParam.epochs       = 500;          % maximum number of epochs to train (for early stopping)
network.trainParam.showWindow   = true;         % show training GUI
network.trainParam.min_grad     = 1e-6;         % minimum performance gradient (for early stopping)
network.trainParam.show         = 1;            % Epochs between displays
network.divideParam.trainRatio  = 0.9;            % all data used for training
network.divideParam.valRatio    = 0.1;            % all data used for training

[net, record] = train(network,X,Y);     % train network

save('net_AEdeep.mat','net')
save('record_AE.mat','record')
