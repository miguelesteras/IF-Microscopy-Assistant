%% Phase 5.4 Train an FeedForward Network from Fourier Descriptors
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a FeedForward network using the fourier descriptors 
%   representation of cell shapes.
%   2. Train an FeedForward Network from Rho descriptor vectors
%   ======================================================================

% 1. Train an FeedForward Network from Fourier Descriptors
clear; close all; clc;
rng('default')

load('fourierDescriptorC4.mat','fourierDescriptor');
Input = fourierDescriptor(:,1:end-2);
Output = fourierDescriptor(:,end-1);

% transform input and output from cell array to column vectors
data = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Input,'UniformOutput',false)))';
output = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Output,'UniformOutput',false)))';

% set apart a random set of samples for validation
valRatio = .15;                        
vIdx = randi([1 size(data,2)],1,round(size(data,2)*valRatio)); 
X = data;    
X(:,vIdx) = [];
Y = output;
Y(:,vIdx) = [];
Xval = data(:,vIdx);
Yval = output(:,vIdx);
target = num2cell(Yval,1);
preFrame = num2cell(Xval(end-15:end,:),1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

% Training feedforward network
hiddenCell = {50, 100, 150, [150 100], [150 100 50], [200 150 50 20]};    
maxEpoch = 1000;
for h = 1:numel(hiddenCell)
    record = zeros(maxEpoch,1);
    hiddenSizes = hiddenCell{h};
    net = feedforwardnet(hiddenSizes);        % create feedforward network
    % Fine-tunning
    net.trainFcn                    = 'trainoss';     % One-step secant backpropagation
    net.divideFcn                   = 'dividerand';   % data division random
    net.divideParam .trainRatio     = 1;              % all data used for training
    net.divideParam.valRatio        = 0;
    net.divideParam.testRatio       = 0;
    net.trainParam.epochs           = 1;           % maximum number of epochs to train (for early stopping)
    net.trainParam.showCommandLine  = false;          % do not generate command-line output
    net.trainParam.showWindow       = true;           % show training GUI
    net.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)

    for e = 1:maxEpoch     
        % train net
        net = train(net,X,Y);        
        % predict
        output = num2cell(net(Xval),1);
        record(e) = sqrt(mse(gsubtract(target,output)));
    end
    record = mat2str([record ones(size(record))*basePerf]);    
    
    save(strcat('net_feedfor_Fourier',num2str(h),'.mat'),'net')
    save(strcat('perf_feedfor_Fourier',num2str(h),'.mat'),'record') 
end


%% 2. Train an FeedForward Network from Rho Descriptors

clear;

load('RhoDescriptors.mat','RhoDescriptors');
Input = RhoDescriptors(:,1:end-1);
Output = RhoDescriptors(:,end);

% transform input and output from cell array to column vectors
data = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Input,'UniformOutput',false)))';
output = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Output,'UniformOutput',false)))';

% set apart a random set of samples for validation
valRatio = .15;                        
vIdx = randi([1 size(data,2)],1,round(size(data,2)*valRatio)); 
X = data;    
X(:,vIdx) = [];
Y = output;
Y(:,vIdx) = [];
Xval = data(:,vIdx);
Yval = output(:,vIdx);
target = num2cell(Yval,1);
preFrame = num2cell(Xval(end-15:end,:),1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

% Training feedforward network
hiddenCell = {50, 100, 150, [150 100], [150 100 50], [200 150 50 20]};    
maxEpoch = 1000;
for h = 1:numel(hiddenCell)
    record = zeros(maxEpoch,1);
    hiddenSizes = hiddenCell{h};
    net = feedforwardnet(hiddenSizes);        % create feedforward network
    % Fine-tunning
    net.trainFcn                    = 'trainoss';     % One-step secant backpropagation
    net.divideFcn                   = 'dividerand';   % data division random
    net.divideParam .trainRatio     = 1;              % all data used for training
    net.divideParam.valRatio        = 0;
    net.divideParam.testRatio       = 0;
    net.trainParam.epochs           = 1;           % maximum number of epochs to train (for early stopping)
    net.trainParam.showCommandLine  = false;          % do not generate command-line output
    net.trainParam.showWindow       = true;           % show training GUI
    net.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)

    for e = 1:maxEpoch     
        % train net
        net = train(net,X,Y);
        % predict
        output = num2cell(net(Xval),1);
        record(e) = sqrt(mse(gsubtract(target,output)));
    end
    record = mat2str([record ones(size(record))*basePerf]);    
    
    save(strcat('net_feedfor_Rho',num2str(h),'.mat'),'net')
    save(strcat('perf_feedfor_Rho',num2str(h),'.mat'),'record') 
end


%% 2. Train an FeedForward Network from Boundary Descriptors

clear;

load('BoundaryDescriptors.mat','BoundaryDescriptors');
Input = BoundaryDescriptors(:,1:end-1);
Output = BoundaryDescriptors(:,end);

% transform input and output from cell array to column vectors
data = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Input,'UniformOutput',false)))';
output = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Output,'UniformOutput',false)))';

% set apart a random set of samples for validation
valRatio = .15;                        
vIdx = randi([1 size(data,2)],1,round(size(data,2)*valRatio)); 
X = data;    
X(:,vIdx) = [];
Y = output;
Y(:,vIdx) = [];
Xval = data(:,vIdx);
Yval = output(:,vIdx);
target = num2cell(Yval,1);
preFrame = num2cell(Xval(end-39:end,:),1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

% Training feedforward network
hiddenCell = {100, 200, 300, [300 200], [300 200 100], [400 300 100 60]};    
maxEpoch = 1000;
for h = 1:numel(hiddenCell)
    record = zeros(maxEpoch,1);
    hiddenSizes = hiddenCell{h};
    net = feedforwardnet(hiddenSizes);        % create feedforward network
    % Fine-tunning
    net.trainFcn                    = 'trainoss';     % One-step secant backpropagation
    net.divideFcn                   = 'dividerand';   % data division random
    net.divideParam .trainRatio     = 1;              % all data used for training
    net.divideParam.valRatio        = 0;
    net.divideParam.testRatio       = 0;
    net.trainParam.epochs           = 1;           % maximum number of epochs to train (for early stopping)
    net.trainParam.showCommandLine  = false;          % do not generate command-line output
    net.trainParam.showWindow       = true;           % show training GUI
    net.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)

    for e = 1:maxEpoch     
        % train net
        net = train(net,X,Y);
        % predict
        output = num2cell(net(Xval),1);
        record(e) = sqrt(mse(gsubtract(target,output)));
    end
    record = mat2str([record ones(size(record))*basePerf]);    
    
    save(strcat('net_feedfor_Boundary',num2str(h),'.mat'),'net')
    save(strcat('perf_feedfor_Boundary',num2str(h),'.mat'),'record') 
end
