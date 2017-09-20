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
clc;clear;close all
rng('default')

load('RhoDescriptors.mat','RhoDescriptors');
data = RhoDescriptors;
seqLen = size(data,2);
samples = size(data,1);
maxSam = floor(samples*0.85);

% define recurrent neural net of type NARNET
trainFcn = 'trainoss';
feedbackDelays = 1:seqLen-1;
hiddenLayerSize = 100;
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

% set apart a random set of samples for validation
valRatio = .15;                        
vIdx = randi([1 samples],1,round(samples*valRatio)); 
trainSet = data;    
trainSet(vIdx,:) = []; 
trainSet = trainSet(randperm(size(trainSet,1)),:);
valSet = data(vIdx,:);
target = valSet(:,end);

% train different model with a dataset of increasing size
epochs = [20 50 100 200 400 round(linspace(600,maxSam,11))];
netPerf  = zeros(numel(epochs),1);

start = 1;
for m = 1:numel(epochs)
    T = trainSet(start:epochs(m),:);
    for i = 1:size(T,1)
        [x,xi,ai,t] = preparets(net,{},{},T(i,:));
        net = train(net,x,t,xi,ai);
    end
    start = epochs(m)+1;
    
    % precit and performance
    output = [];
    for j = 1:size(valSet,1)
        query = valSet(j,1:end-1);
        [xq,xiq,aiq,tq] = preparets(net,{},{},query);
        output = [output ; net(xq,xiq,aiq)];
    end
    netPerf(m) = sqrt(mse(gsubtract(target,output)));
end

preFrame = valSet(:,end-1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

rho_performance = mat2str([netPerf ones(size(netPerf))*basePerf]);
save('perf_Rho.mat','rho_performance')

%% Increasing number of epochs

epochs = 50;
netPerf  = zeros(epochs,1);
for m = 1:50
    for i = 1:size(trainSet,1)
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
    netPerf(m) = sqrt(mse(gsubtract(target,output)));
end

rho_performance = mat2str([netPerf ones(size(netPerf))*basePerf]);
save('perf_Rho_1to50Epochs.mat','rho_performance')
save('net_Rho_narnet_50epochs.mat','net')




% %% old 
% for m = 1:numel(epochs)
%     
%     % set apart a random set of samples for validation
%     valRatio = .15;                        
%     vIdx = randi([1 samples],1,round(samples*valRatio)); 
%     trainSet = data;    
%     trainSet(vIdx,:) = []; 
%     trainSet = trainSet(randperm(size(trainSet,1)),:);
%     valSet   = data(vIdx,:);
%     
%     net = narnet(feedbackDelays,hiddenLayerSize,'open',trainFcn);
%     net = removedelay(net);
%     net.trainParam.epochs   = 2;
%     net.divideFcn           = 'dividerand'; 
%     net.divideMode          = 'time';
%     net.performFcn          = 'mse';  
%     net.input.processFcns   = {'removeconstantrows','mapminmax'};
%     net.divideParam.trainRatio  = 1;
%     net.trainParam.showWindow   = false;           
%     net.trainParam.min_grad = 0;
%     
%     for i = 1:epochs(m)
%         T = trainSet(i,:);
%         [x,xi,ai,t] = preparets(net,{},{},T);
%         net = train(net,x,t,xi,ai);
%     end
%     
%     % precit and performance
%     output = [];
%     for j = 1:size(valSet,1)
%         query = valSet(j,1:end-1);
%         [xq,xiq,aiq,tq] = preparets(net,{},{},query);
%         output = [output ; net(xq,xiq,aiq)];
%     end
%     target = valSet(:,end);
%     preFrame = valSet(:,end-1);
%     netPerfor(m) = sqrt(mse(gsubtract(target,output)));
%     basePerfor(m) = sqrt(mse(gsubtract(target,preFrame)));
% end
% rho_performance = mat2str([netPerfor basePerfor]);
% save('perf_Rho_original.mat','rho_performance')
% 
% %% Increasing epochs
% 
% epochs = linspace(1,10,10);
% netPerfor  = zeros(numel(epochs),1);
% basePerfor = zeros(numel(epochs),1);
% 
% for m = 1:numel(epochs)
%     
%     % set apart a random set of samples for validation
%     valRatio = .15;                        
%     vIdx = randi([1 samples],1,round(samples*valRatio)); 
%     trainSet = data;    
%     trainSet(vIdx,:) = []; 
%     trainSet = trainSet(randperm(size(trainSet,1)),:);
%     valSet = data(vIdx,:);
%     
%     net = narnet(feedbackDelays,hiddenLayerSize,'open',trainFcn);
%     net = removedelay(net);
%     net.trainParam.epochs   = 2;
%     net.divideFcn           = 'dividerand'; 
%     net.divideMode          = 'time';
%     net.performFcn          = 'mse';  
%     net.input.processFcns   = {'removeconstantrows','mapminmax'};
%     net.divideParam.trainRatio  = 1;
%     net.trainParam.showWindow   = false;           
%     net.trainParam.min_grad = 0;
%     % number of epochs to train with full training set, from 1 to 10.
%     for k = 1:epochs(m)
%         for i = 1:size(trainSet,1)
%             T = trainSet(i,:);
%             [x,xi,ai,t] = preparets(net,{},{},T);
%             net = train(net,x,t,xi,ai);
%         end
%     end    
%     % precit and performance
%     output = [];
%     for j = 1:size(valSet,1)
%         query = valSet(j,1:end-1);
%         [xq,xiq,aiq,tq] = preparets(net,{},{},query);
%         output = [output ; net(xq,xiq,aiq)];
%     end
%     target = valSet(:,end);
%     preFrame = valSet(:,end-1);
%     netPerfor(m) = sqrt(mse(gsubtract(target,output)));
%     basePerfor(m) = sqrt(mse(gsubtract(target,preFrame)));
% end
% rho_performance = mat2str([netPerfor basePerfor]);
% save('perf_Rho_epochs.mat','rho_performance')
