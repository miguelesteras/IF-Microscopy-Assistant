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

seqLen = size(Input,2);
samples = size(Input,1);
maxSam = floor(samples*0.85);

% transform input and output from cell array to column vectors
X = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Input,'UniformOutput',false)))';
Y = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Output,'UniformOutput',false)))';

% Training feedforward network
hiddenCell = {100, [150 100], [150 100 50]};    

for h = 1:numel(hiddenCell)
    hiddenSizes = hiddenCell{h};
    net = feedforwardnet(hiddenSizes);        % create feedforward network
    % Fine-tunning
    net.trainFcn                    = 'trainoss';     % One-step secant backpropagation
    net.divideFcn                   = 'dividerand';   % data division random
    net.divideParam .trainRatio     = 0.85;              % all data used for training
    net.divideParam.valRatio        = 0.15;
    net.divideParam.testRatio       = 0;
    net.trainParam.epochs           = 100;           % maximum number of epochs to train (for early stopping)
    net.trainParam.showCommandLine  = false;          % do not generate command-line output
    net.trainParam.showWindow       = true;           % show training GUI
    net.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)
    net.trainParam.show             = 1;              % Epochs between displays
    % train net
    [net,record] = train(net,X,Y);
    % save network and performance
    save(strcat('net_feedfor_Fourier',num2str(h),'.mat'),'net')
    save(strcat('perf_feedfor_Fourier',num2str(h),'.mat'),'record')
end


% 2. Train an FeedForward Network from Rho Descriptors

clear;

load('RhoDescriptors.mat','RhoDescriptors');
Input = RhoDescriptors(:,1:end-1);
Output = RhoDescriptors(:,end);

seqLen = size(Input,2);
samples = size(Input,1);
maxSam = floor(samples*0.85);

% transform input and output from cell array to column vectors
X = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Input,'UniformOutput',false)))';
Y = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),Output,'UniformOutput',false)))';

% Training feedforward network
hiddenCell = {100, [150 100], [150 100 50]};    

for h = 1:numel(hiddenCell)
    hiddenSizes = hiddenCell{h};
    net = feedforwardnet(hiddenSizes);        % create feedforward network
    % Fine-tunning
    net.trainFcn                    = 'trainoss';     % One-step secant backpropagation
    net.divideFcn                   = 'dividerand';   % data division random
    net.divideParam .trainRatio     = 0.85;              % all data used for training
    net.divideParam.valRatio        = 0.15;
    net.divideParam.testRatio       = 0;
    net.trainParam.epochs           = 100;           % maximum number of epochs to train (for early stopping)
    net.trainParam.showCommandLine  = false;          % do not generate command-line output
    net.trainParam.showWindow       = true;           % show training GUI
    net.trainParam.max_fail         = 10;             % maximum validation failures (for early stopping)
    net.trainParam.show             = 1;              % Epochs between displays
    % train net
    [net,record] = train(net,X,Y);
    % save network and performance
    save(strcat('net_feedfor_Rho',num2str(h),'.mat'),'net')
    save(strcat('perf_feedfor_Rho',num2str(h),'.mat'),'record')
end