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

load('fourierDescriptorC4.mat','fourierDescriptor');
data = fourierDescriptor;

seqLen = size(data,2);
samples = size(data,1);

trainFcn = 'trainoss';
feedbackDelays = 1:seqLen-2;
hiddenLayerSize = 100;


epochs = [10 20 50 100 500 linspace(1000, 15000, 29)];
netPerfor  = zeros(numel(epochs),1);
basePerfor = zeros(numel(epochs),1);
for m = 1:numel(epochs)
    [tIdx, vIdx] = dividerand(samples, .85, .15);
    trainSet = [data(tIdx,:) ; data(tIdx,:) ;
        data(tIdx,:) ; data(tIdx,:) ; data(tIdx,:) ; data(tIdx,:);
        data(tIdx,:) ; data(tIdx,:) ; data(tIdx,:) ; data(tIdx,:)];        
    valSet   = data(vIdx,:);
    for i = 1:epochs(m)
        T = trainSet(i,:);
        
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
    netPerfor(m)  = perform(net,target,output)/numel(target);
    basePerfor(m) = perform(net,target,preFrame)/numel(target);
end

fourier_performance = [netPerfor basePerfor];

save('net_fourier_narnet.mat','net')
save('perf_fourier.mat','fourier_performance')
