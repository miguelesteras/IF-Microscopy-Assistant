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

rng('default')

%% Rho

load('RhoDescriptors.mat','RhoDescriptors');
load('RhoDescriptors04.mat','RhoDescriptors04');
load('RhoDescriptors05.mat','RhoDescriptors05');
load('RhoDescriptors06.mat','RhoDescriptors06');
load('RhoDescriptors07.mat','RhoDescriptors07');
data = RhoDescriptors;
data_syn = [RhoDescriptors04 ; RhoDescriptors05 ;
        RhoDescriptors06 ; RhoDescriptors07];

seqLen = size(data,2);
samples = size(data,1);
maxSam = floor(samples*0.85)+size(data_syn,1);

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
trainSet = [trainSet ; data_syn];
trainSet = trainSet(randperm(size(trainSet,1)),:);
valSet = data(vIdx,:);
target = valSet(:,end);
preFrame = valSet(:,end-1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

epochs = 100;
netPerf  = zeros(epochs,1);
for m = 1:epochs
    % check progress
    save('progress_narnet100_rho.mat','m','-ascii')

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

record = mat2str([netPerf ones(size(netPerf))*basePerf]);
save('perf_Rho_Narnet100.mat','record')
save('net_Rho_Narnet100.mat','net')

%% Fourier

load('fourierDescriptorC4.mat','fourierDescriptor');
data = fourierDescriptor(:,1:end-1);

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
preFrame = valSet(:,end-1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

epochs = 100;
netPerf  = zeros(epochs,1);
for m = 1:epochs
    % check progress
    save('progress_narnet100_fourier.mat','m','-ascii')

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

record = mat2str([netPerf ones(size(netPerf))*basePerf]);
save('perf_Fourier_Narnet100.mat','record')
save('net_Fourier_Narnet100.mat','net')

%% Boundary

load('BoundaryDescriptors.mat','BoundaryDescriptors');
data = cellfun(@(x) (reshape(x, [], 1)),BoundaryDescriptors,'UniformOutput',false);

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
preFrame = valSet(:,end-1);
basePerf = sqrt(mse(gsubtract(target,preFrame)));

epochs = 100;
netPerf  = zeros(epochs,1);
for m = 1:epochs
    % check progress
    save('progress_narnet100_boundary.mat','m','-ascii')

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

record = mat2str([netPerf ones(size(netPerf))*basePerf]);
save('perf_Boundary_Narnet100.mat','record')
save('net_Boundary_Narnet100.mat','net')