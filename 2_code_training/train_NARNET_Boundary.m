%% Phase 5.7 Train Recurrent Network NARNET (Matlab)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a NARNET recurrent network using Fourier descriptors, boundary
%   descriptors and Rho descriptors.
%   ======================================================================
clear; close all; clc;
rng('default')

load('BoundaryDescriptors.mat','BoundaryDescriptors');
data = cellfun(@(x) (reshape(x, [], 1)),BoundaryDescriptors,'UniformOutput',false);

seqLen = size(data,2);
samples = size(data,1);

trainFcn = 'trainoss';
feedbackDelays = 1:seqLen-1;
hiddenLayerSize = 200;


epochs = [10 20 50 100 200 400 linspace(600,floor(samples*0.85),10)];
netPerfor  = zeros(numel(epochs),1);
basePerfor = zeros(numel(epochs),1);
for m = 1:numel(epochs)
    valRatio = .15;                        
    vIdx = randi([1 samples],1,round(samples*valRatio)); 
    trainSet = data;    
    trainSet(vIdx,:) = []; 
    trainSet = trainSet(randperm(size(trainSet,1)),:);
    valSet   = data(vIdx,:);
    
    net = narnet(feedbackDelays,hiddenLayerSize,'open',trainFcn);
    net = removedelay(net);
    net.trainParam.epochs   = 2;
    net.divideFcn           = 'dividerand'; 
    net.divideMode          = 'time';
    net.performFcn          = 'mse';  
    net.input.processFcns   = {'removeconstantrows','mapminmax'};
    net.divideParam.trainRatio  = 1;
    net.trainParam.showWindow   = false;           
    net.trainParam.min_grad = 0;
           
    for i = 1:epochs(m)
        T = trainSet(i,:);
        [x,xi,ai,t] = preparets(net,{},{},T);
        net = train(net,x,t,xi,ai);
    end
    
    % precit and performance
    output = [];
    for j = 1:size(valSet,1)
        query = valSet(j,1:end-1);
        [xq,xiq,aiq,tq] = preparets(net,{},{},query);
        output = [output ; net(xq,xiq,aiq)];
    end
    target = valSet(:,end);
    preFrame = valSet(:,end-1);
    netPerfor(m) = sqrt(mse(gsubtract(target,output)));
    basePerfor(m) = sqrt(mse(gsubtract(target,preFrame)));
end

boundary_performance = mat2str([netPerfor basePerfor]);
save('perf_Boundary.mat','boundary_performance')

%% Increasing number of Epochs

epochs = linspace(1,10,10);
netPerfor  = zeros(numel(epochs),1);
basePerfor = zeros(numel(epochs),1);
for m = 1:numel(epochs)
    valRatio = .15;                        
    vIdx = randi([1 samples],1,round(samples*valRatio)); 
    trainSet = data;    
    trainSet(vIdx,:) = []; 
    trainSet = trainSet(randperm(size(trainSet,1)),:);
    valSet   = data(vIdx,:);
    
    net = narnet(feedbackDelays,hiddenLayerSize,'open',trainFcn);
    net = removedelay(net);
    net.trainParam.epochs   = 2;
    net.divideFcn           = 'dividerand'; 
    net.divideMode          = 'time';
    net.performFcn          = 'mse';  
    net.input.processFcns   = {'removeconstantrows','mapminmax'};
    net.divideParam.trainRatio  = 1;
    net.trainParam.showWindow   = false;           
    net.trainParam.min_grad = 0;
    % number of epochs to train with full training set, from 1 to 10.
    for k = 1:epochs(m)       
        for i = 1:size(trainSet,1)
            T = trainSet(i,:);
            [x,xi,ai,t] = preparets(net,{},{},T);
            net = train(net,x,t,xi,ai);
        end
    end
    
    % precit and performance
    output = [];
    for j = 1:size(valSet,1)
        query = valSet(j,1:end-1);
        [xq,xiq,aiq,tq] = preparets(net,{},{},query);
        output = [output ; net(xq,xiq,aiq)];
    end
    target = valSet(:,end);
    preFrame = valSet(:,end-1);
    netPerfor(m) = sqrt(mse(gsubtract(target,output)));
    basePerfor(m) = sqrt(mse(gsubtract(target,preFrame)));
end

boundary_performance = mat2str([netPerfor basePerfor]);
save('net_Boundary_narnet.mat','net')
save('perf_Boundary_epochs.mat','boundary_performance')
