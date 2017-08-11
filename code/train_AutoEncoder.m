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

rng('default')

% Create a joined dataSet from all video files
files = dir('*_metadata.mat');      
num_files = length(files);
inputData = [];   
targetData = [];
for i = 1:num_files
    load(files(i).name,'metadata');                               
    load(strcat(metadata.name,'_dynamicInput.mat'),'dynamicInput');     
    load(strcat(metadata.name,'_dynamicTarget.mat'),'dynamicTarget');     
    inputData = [inputData dynamicInput'];
    targetData = [targetData dynamicTarget'];
end

%% reduce size of images to reduce size of Autoencoder
load(strcat(metadata.name,'_maxRadii.mat'),'maxRadii');     
index = find(~cellfun(@isempty,maxRadii)); % non empty cells in cell sequences array
radious = ceil(cellfun(@(v) v(2), maxRadii(index)));    
%histogram(radious)

cropInput = cellfun(@(x) imcrop(x,[50 35 93 123]),inputData,...
    'UniformOutput',false);

reduceTarget = cellfun(@(x) imresize(imcrop(x,[50 35 93 123]),0.5),targetData,...
    'UniformOutput',false);

%% Train a autoencoder for dynamic images: First layer
enco_fnc = ['logsig' 'satlin'];
deco_fnc = ['logsig' 'satlin' 'purelin'];

hiddenSize1 = 1000; tic;                         	% start timer
AE1 = trainAutoencoder(cropInput,hiddenSize1, ... % images in cell array format
    'MaxEpochs',1000, ...                       % max no.epochs for early stopping
    'EncoderTransferFunction', 'logsig',...     % logistic sigmoid function
    'DecoderTransferFunction', 'logsig',...     % logistic sigmoid function
    'LossFunction', 'msesparse',...             % loss function used for training
    'TrainingAlgorithm', 'trainscg',...         % training algorithm, scaled conjugate gradient descent
    'L2WeightRegularization',0.004, ...         % coefficient for the L2 weight regularizer in the cost function
    'SparsityRegularization',4, ...             % controls the weighting of the sparsity regularizer
    'SparsityProportion',0.2, ...               % proportion of training examples which a neuron in the hidden layer should activate in response to.
    'ScaleData', false);                        % do not rescale input data

view(AE1)                                  % visualize autoencoder1
plotWeights(AE1);                          % visualize feautres from autoencoder1 hidden layer
sprintf('Training time AE1: %0.0f %s',toc,'s')  % duration of training

% i = 100;
% Reconstructed = predict(autoenc1,Xtrain{i});    % code -> decode image 1 to 10
% figure; subplot(1,2,1); imshow(Xtrain{i});      % show original image
% subplot(1,2,2); imshow(Reconstructed);          % show reconstracted image

AE1features = encode(AE1,cropInput);            % generate feautres from autoencoder1 hidden layer

%% Second layer

% use features extracted by autoencoder1 as input to autoencoder 2

hiddenSize2 = 500; tic;                                  % start timer
AE2 = trainAutoencoder(AE1features,hiddenSize2, ... 
    'MaxEpochs',1000, ...
    'EncoderTransferFunction', 'logsig',...
    'DecoderTransferFunction', 'logsig',...
    'LossFunction', 'msesparse',...
    'TrainingAlgorithm', 'trainscg',...                 % training algorithm, scaled conjugate gradient descent
    'L2WeightRegularization',0.004, ...                 % coefficient for the L2 weight regularizer in the cost function
    'SparsityRegularization',4, ...
    'SparsityProportion',0.2, ...
    'ScaleData', false);                                % do not rescale input data

sprintf('Training time AutoEncoder2: %0.0f %s',toc,'s') % duration of training
view(AE2)
AE2features = encode(AE2,AE1features);

%% Train a fully connected layer

% transform input and output from cell array to column vectors
X = cell2mat(cellfun(@(x) (reshape(x, 1, []))',cropInput,'UniformOutput',false));
Y = cell2mat(cellfun(@(x) (reshape(x, 1, []))',reduceTarget,'UniformOutput',false));

% train network stack
FC = fullyConnectedLayer(2914);
network = stack(AE1,AE2,FC);
view(network)

network.inputs{1}.size = 11656;
network.trainFcn                = 'trainoss';   % One-step secant backpropagation
network.trainParam.epochs       = 100;          % maximum number of epochs to train (for early stopping)
network.trainParam.showWindow   = true;         % show training GUI
network.trainParam.min_grad     = 1e-6;         % minimum performance gradient (for early stopping)
network.trainParam.show         = 1;            % Epochs between displays
network.divideParam.trainRatio = 0.9;            % all data used for training
network.divideParam.valRatio = 0.1;            % all data used for training

[trainNet, record] = train(network,X,Y);     % train network



