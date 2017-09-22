%% Extras. Create Plots for feedforward network
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Boundary
clear; close all; clc

Pfiles = dir('perf_feedfor_Boundary*'); 
Nfiles = dir('net_feedfor_Boundary*');
num_files = length(Pfiles);
Perf = cell(num_files,2);
Nets = cell(num_files,2);

for i = 1:num_files
    load(Pfiles(i).name,'record');                                     
    Perf{i,1} = str2num(record);
    Perf{i,2} = Pfiles(i).name;
    
    load(Nfiles(i).name,'net');
    Nets{i,1} = net;
    Nets{i,2} = Nfiles(i).name; 
end
  
% Plot performance (validation error during training)
epochs = size(Perf{1,1},1);
x = linspace(1,epochs,epochs);
base = Perf{1,1}(:,2);

figure
plot(x, base, 'Color', [0 0 .5], 'LineWidth', 2); hold on
for j = 1:num_files
    y = Perf{j,1}(:,1);
    plot(x, y, 'LineWidth', 3);
end

xlim([min(x) max(x)])
grid on
xlabel('Epochs')
ylabel('Root Mean Square Error (log scale)')
yticks([10 linspace(20,50,4) 100])
set(gca,'fontsize',16,'YScale', 'log')
lg = legend('Baseline','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6');
lg.FontSize = 16;

% Show reconstruction of (numEx) examples 
numEx = 2;
load('BoundaryDescriptors.mat','BoundaryDescriptors');
dataSet = BoundaryDescriptors;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
% convert cell arrays into matrices (columns = samples)
input = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,1:end-1),'UniformOutput',false)))';
target = dataSet(idx,end); 
preFrame = dataSet(idx,end-1);

net = Nets{5,1};
prediction = net(input);

load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);

for k = 1:numel(idx)    
    selection = cells{k};
    selection(:,1) = wrapTo2Pi(selection(:,1) + pi); 
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
    
    mask2 = (target{k}*diag([2 -2]))+100;
    I2 = insertMarker(I,mask2,'Size',5, 'Color','red');
    figure
    imshow(I2)
    mask3 = (preFrame{k}*diag([2 -2]))+100;
    I3 = insertMarker(I,mask3,'Size',5, 'Color','blue');
    figure
    imshow(I3)
    i1 = linspace(1,39,20);
    i2 = linspace(2,40,20);
    mask4 = [prediction(i1,k) prediction(i2,k)];
    mask4 = (mask4*diag([2 -2]))+100;
    I4 = insertMarker(I,mask4,'Size',5, 'Color','blue');
    figure
    imshow(I4)    
end

%% Dice coefficient

numEx = 2;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);

% convert cell arrays into matrices (columns = samples)
input = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,1:end-1),'UniformOutput',false)))';
target = dataSet(idx,end); 
preFrame = dataSet(idx,end-1);

net = Nets{5,1};
prediction = num2cell(net(input),1);
% reshape prediction column matrix into cell array 
i1 = linspace(1,39,20);
i2 = linspace(2,40,20);
prediction = (cellfun(@(x) ([x(i1) x(i2)]),prediction,'UniformOutput',false));


load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);
dice = zeros(3,numel(idx));

for k = 1:numel(idx)    
    selection = cells{k};
    selection(:,1) = wrapTo2Pi(selection(:,1) + pi); 
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
        points = (vectors{v}{k}*diag([2 -2]))+100;
        canvas2 = poly2mask(points(:,1),points(:,2),200,200);
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