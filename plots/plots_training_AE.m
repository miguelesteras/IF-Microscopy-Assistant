%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

clear; close all; clc;
rng('default')

load('dynamicImages.mat','dynamicImages');
load('dynamicImages04.mat','dynamicImages04');
load('dynamicImages05.mat','dynamicImages05');
load('dynamicImages06.mat','dynamicImages06');
load('dynamicImages07.mat','dynamicImages07');

dynamicImagesAll = [dynamicImages ; dynamicImages04 ; dynamicImages05 ;
            dynamicImages06 ; dynamicImages07];
        
inputData = cell(size(dynamicImagesAll));
% Compute a Gaussian pyramid reduction by one level
for j = 1:size(dynamicImagesAll,1)
    inputData{j,1} = impyramid(dynamicImagesAll{j,1}, 'reduce');    
end

dynamicMini = cell(size(dynamicImages));
% Compute a Gaussian pyramid reduction by one level
for j = 1:size(dynamicImages,1)
    dynamicMini{j,1} = impyramid(dynamicImages{j,1}, 'reduce');    
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

%% AE1 performance
load('/Users/miguelesteras/Desktop/matlab_project/AE/net_AE1.mat','AE1');                                     
%view(AE1)                                  % visualize autoencoder1
%plotWeights(AE1);                          % visualize feautres from autoencoder1 hidden layer

for i = 300:302
    Reconstructed = predict(AE1,Xtrain{i});         % code -> decode image 1 to 10
    figure; subplot(1,2,1); imshow(Xtrain{i});      % show original image
    subplot(1,2,2); imshow(Reconstructed);          % show reconstracted image
end
%% AE2 performance
load('/Users/miguelesteras/Desktop/matlab_project/AE/net_AE2.mat','AE2');                                     
%view(AE2)                                  
%plotWeights(AE2);                       

%% feedforward layer
load('/Users/miguelesteras/Desktop/matlab_project/AE/net_AEfeedfor.mat','net');                                     
load('/Users/miguelesteras/Desktop/matlab_project/AE/perf_AEfeedfor.mat','record');                                     
%view(net)                                  

% Plot performance (validation error during training)
x = record.epoch;
y = record.perf;
[y1, x1] = min(record.perf);

figure
plot(x, y, 'LineWidth', 3); hold on
scatter(x1, y1,100,'o','r','MarkerFaceAlpha',0.7,'MarkerFaceColor',[1 0.3 0.3]);
ylim([0 1000])
xlim([min(x) max(x)+5])
grid on
xlabel('Epochs')
ylabel('Mean Square Error')
set(gca,'fontsize',16, 'YScale','log')
lg = legend('Validation Performance',strcat('min. value = ',num2str(min(record.perf))));
lg.FontSize = 16;

%% Deep network performance

load('/Users/miguelesteras/Desktop/matlab_project/AE/net_AEdeep.mat','net');                                     
load('/Users/miguelesteras/Desktop/matlab_project/AE/record_AE.mat','record');                                     
%view(net)                                  

% Plot performance (validation error during training)
x = record.epoch;
y = record.perf;
[y1, x1] = min(record.perf);

figure
plot(x, y, 'LineWidth', 3, 'Color', [0 .5 0]); hold on
scatter(x1, y1,100,'o','r','MarkerFaceAlpha',0.7,'MarkerFaceColor',[1 0.3 0.3]);
ylim([0 1000])
xlim([min(x) max(x)+5])
grid on
xlabel('Epochs')
ylabel('Mean Square Error')
set(gca,'fontsize',16, 'YScale','log')
lg = legend('Validation Performance',strcat('min. value = ',num2str(min(record.perf))));
lg.FontSize = 16;


%% Show reconstruction of (numEx) examples  

load('RhoDescriptors.mat','RhoDescriptors');
load('dynamicImages.mat','dynamicImages');
load('/Users/miguelesteras/Desktop/matlab_project/AE/net_AEdeep.mat','net');                                     

dynamicMini = cell(size(dynamicImages,1),1);
% Compute a Gaussian pyramid reduction by one level
for j = 1:size(dynamicImages,1)
    dynamicMini{j,1} = impyramid(dynamicImages{j,1}, 'reduce');    
end

numEx = 2;
dataSet = RhoDescriptors;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);

% convert cell arrays into matrices (columns = samples)
input = cell2mat(cellfun(@(x) (reshape(x, [], 1)),dynamicMini(idx,1)','UniformOutput',false));
target = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end),'UniformOutput',false)))';
preFrame = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end-1),'UniformOutput',false)))';

prediction = net(input);

vecSize   = size(prediction,1);            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);

for k = 1:numel(idx)    
    figure;
    imshow(dynamicImages{idx(k),1});
    figure;
    polarscatter(cells{k}(:,1), cells{k}(:,2),'filled',...
        'MarkerFaceAlpha',0.7,'MarkerFaceColor',[.7 .7 .7],...
        'MarkerEdgeColor','none'); hold on;
    polarscatter(refAngles,target(:,k), 'SizeData',100);
    for m = 1:numel(target(:,k))
        polarplot([0;refAngles(m)],[0;target(m,k)],'Color','red');
    end
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    set(gca,'fontsize',16)
    rticks([])
    legend('original shape','Rho features')
    hold off;
    
    figure;
    polarscatter(cells{k}(:,1), cells{k}(:,2),...
        'MarkerFaceAlpha',0.7,'MarkerFaceColor',[.7 .7 .7],...
        'MarkerEdgeColor','none'); hold on;
    rhoLine = [refAngles' preFrame(:,k)];
    rhoLine = [rhoLine ; rhoLine(1,:)];
    polarplot(rhoLine(:,1),rhoLine(:,2),'b','lineWidth',3); 
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    set(gca,'fontsize',16)
    rticks([])
    legend('original shape','previous frame')
    hold off

    figure;
    polarscatter(cells{k}(:,1), cells{k}(:,2),...
        'MarkerFaceAlpha',0.7,'MarkerFaceColor',[.7 .7 .7],...
        'MarkerEdgeColor','none'); hold on;
    rhoLine = [refAngles' prediction(:,k)];
    rhoLine = [rhoLine ; rhoLine(1,:)];
    polarplot(rhoLine(:,1),rhoLine(:,2),'g','lineWidth',3); 
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    set(gca,'fontsize',16)
    rticks([])
    legend('original shape','model prediction')
    hold off
end

%% Dice coefficient

numEx = 200;
dataSet = RhoDescriptors;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);

% convert cell arrays into matrices (columns = samples)
input = cell2mat(cellfun(@(x) (reshape(x, [], 1)),dynamicMini(idx,1)','UniformOutput',false));
target = data(:,end);
preFrame = data(:,end-1);

prediction = num2cell(net(input),1);

load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);
dice = zeros(3,numel(idx));

for k = 1:numel(idx)    
    selection = cells{k};
    [x,y] = pol2cart(selection(:,1),selection(:,2));
    mask = [round(x) round(y)];
    ind = sub2ind([300,300], (mask(:,2)*-2)+150, (mask(:,1)*2)+150);
    canvas = false(300,300);
    canvas(ind) = true;    
    canvas = imfill(bwmorph(imfill(canvas,'holes'),'bridge'),'holes');
    se = strel('disk',7);
    canvas = imopen(canvas,se);
%     I = double(cat(3, canvas, canvas, canvas));
%     I(I==1) = 0.75;
%     I(I==0) = 1;
    ImgSize = size(canvas);
    vectors = {target preFrame prediction'}; 
    for v = 1:numel(vectors)
        [x,y] = pol2cart(refAngles',vectors{v}{k});
        canvas2 = poly2mask((round(x)*2)+150,(round(y)*-2)+150,300,300);
        red = zeros(ImgSize);
        red(bwmorph(canvas2,'remove')==1) = 255;
%         I2 = I;
%         I2(:,:,1) = I(:,:,1) + red;
%         I2(:,:,2) = I(:,:,2) - red;
%         I2(:,:,3) = I(:,:,3) - red;
%         figure
%         imshow(I2);
        % compute Sørensen-Dice Coefficient between original image and reconstructed 
        dice(v,k) = 2*nnz(canvas&canvas2)/(nnz(canvas2) + nnz(canvas));  
    end
end
save('dice_AE.mat','dice')
