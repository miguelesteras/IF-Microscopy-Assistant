%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Rho descriptors 

clc; clear;
cd '/Users/miguelesteras/Desktop/matlab_project/Narnet/'

load('perf_Rho_Narnet.mat','record');                                     
Perf = str2num(record);

% Plot performance (validation error during training)
epochs = size(Perf,1);
x = linspace(1,epochs,epochs);
base = Perf(:,2);

figure
plot(x, base, 'Color', [0 0 .5], 'LineWidth', 2); hold on
y = Perf(:,1);
plot(x, y, 'LineWidth', 3);
ylim([0 10])
xlim([min(x) max(x)])
grid on
xlabel('Epochs')
ylabel('Root Mean Square Error')
yticks([1 linspace(2,10,5)])
set(gca,'fontsize',16)
lg = legend('Baseline','Model Performance');
lg.FontSize = 16;

% Show reconstruction of (numEx) examples  
numEx = 2;
load('RhoDescriptors.mat','RhoDescriptors');
dataSet = RhoDescriptors;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
target = dataSet(idx,end);
preFrame = dataSet(idx,end-1);
prediction = [];

% precit and performance
load('net_Rho_Narnet.mat','net');
for j = 1:size(data,1)
    query = data(j,1:end-1);
    [xq,xiq,aiq,tq] = preparets(net,{},{},query);
    prediction = [prediction ; net(xq,xiq,aiq)];
end

vecSize   = 16;            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);

for k = 1:numel(idx)    
    figure;
    polarscatter(cells{k}(:,1), cells{k}(:,2),'filled',...
        'MarkerFaceAlpha',0.7,'MarkerFaceColor',[.7 .7 .7],...
        'MarkerEdgeColor','none'); hold on;
    polarscatter(refAngles,target{k}, 'SizeData',100);
    for m = 1:numel(target{k})
        polarplot([0;refAngles(m)],[0;target{k}(m)],'Color','red');
    end
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    set(gca,'fontsize',16)
    rticks([])
    legend('original shape')
    hold off;
    
    figure;
    polarscatter(cells{k}(:,1), cells{k}(:,2),...
        'MarkerFaceAlpha',0.7,'MarkerFaceColor',[.7 .7 .7],...
        'MarkerEdgeColor','none'); hold on;
    rhoLine = [refAngles' preFrame{k}];
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
    rhoLine = [refAngles' prediction{k}];
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
numEx = 2;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
target = dataSet(idx,end);
preFrame = dataSet(idx,end-1);
prediction = [];

% precit and performance
load('net_Rho_Narnet.mat','net');
for j = 1:size(data,1)
    query = data(j,1:end-1);
    [xq,xiq,aiq,tq] = preparets(net,{},{},query);
    prediction = [prediction ; net(xq,xiq,aiq)];
end

load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);
dice = zeros(3,numel(idx));

for k = 1:numel(idx)    
    selection = cells{k};
    [x,y] = pol2cart(selection(:,1),selection(:,2));
    mask = [round(x) round(y)];
    ind = sub2ind([200,200], (mask(:,2)*-2)+100, (mask(:,1)*2)+100);
    canvas = false(200,200);
    canvas(ind) = true;    
    canvas = imfill(bwmorph(imfill(canvas,'holes'),'bridge'),'holes');
    se = strel('disk',7);
    canvas = imopen(canvas,se);
    I = double(cat(3, canvas, canvas, canvas));
    I(I==1) = 0.75;
    I(I==0) = 1;
    ImgSize = size(canvas);
    vectors = {target preFrame prediction}; 
    for v = 1:numel(vectors)
        [x,y] = pol2cart(refAngles',vectors{v}{k});
        canvas2 = poly2mask((round(x)*2)+100,(round(y)*-2)+100,200,200);
        red = zeros(ImgSize);
        red(bwmorph(canvas2,'remove')==1) = 255;
        I2 = I;
        I2(:,:,1) = I(:,:,1) + red;
        I2(:,:,2) = I(:,:,2) - red;
        I2(:,:,3) = I(:,:,3) - red;
        figure
        imshow(I2);
        % compute Sørensen-Dice Coefficient between original image and reconstructed 
        dice(v,k) = 2*nnz(canvas&canvas2)/(nnz(canvas2) + nnz(canvas));  
    end
end
save('dice_narnet_rho.mat','dice')
