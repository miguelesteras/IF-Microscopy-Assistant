%% Extras. Create Plots for feedforward network
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Rho descriptors 

clc; clear;
cd '/Users/miguelesteras/Desktop/matlab_project/feedforward/'
Pfiles = dir('perf_feedfor_Rho*'); 
Nfiles = dir('net_feedfor_Rho*');
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
yticks([1 linspace(2,10,5) 100])
set(gca,'fontsize',16,'YScale', 'log')
lg = legend('Baseline','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6');
lg.FontSize = 16;


% Show reconstruction of (numEx) examples  
numEx = 2;
load('RhoDescriptors.mat','RhoDescriptors');
dataSet = RhoDescriptors;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
% convert cell arrays into matrices (columns = samples)
input = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,1:end-1),'UniformOutput',false)))';
target = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end),'UniformOutput',false)))';
preFrame = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end-1),'UniformOutput',false)))';

net = Nets{5,1};
prediction = net(input);

vecSize   = size(prediction,1);            
refAngles = linspace(0,2*pi,vecSize+1);
refAngles = refAngles(1:end-1);
load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);

for k = 1:numel(idx)    
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
    legend('original shape')
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

%% Boundary

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
numEx = 4;
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


%% Fourier descriptors 

clc; clear;
cd '/Users/miguelesteras/Desktop/matlab_project/feedforward/'
Pfiles = dir('perf_feedfor_Fourier*'); 
Nfiles = dir('net_feedfor_Fourier*');
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
ylim([0 1])
grid on
xlabel('Epochs')
ylabel('Root Mean Square Error')
%yticks([1 linspace(2,10,5) 100])
set(gca,'fontsize',16)
lg = legend('Baseline','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6');
lg.FontSize = 16;


% Show reconstruction of (numEx) examples  
numEx = 2;
load('fourierDescriptorC4.mat','fourierDescriptor');
dataSet = fourierDescriptor;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
% convert cell arrays into matrices (columns = samples)
input = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,1:end-2),'UniformOutput',false)))';
target = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end-1),'UniformOutput',false)))';
preFrame = (cell2mat(cellfun(@(x) (reshape(x', 1, [])),data(:,end-2),'UniformOutput',false)))';

net = Nets{5,1};
prediction = net(input);

load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);
dice = zeros(3,numel(idx));

for k = 1:numel(idx)    
    selection = cells{k};
    selection(:,1) = wrapTo2Pi(selection(:,1) + pi/2); 
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
    
    % inverse fourier transformation
    T = 300;
    s = 0.5;
    ImgSize = size(canvas);
    vectors = {target preFrame prediction}; 
    for v = 1:numel(vectors)        
        a = vectors{v}(1:4,k);
        b = vectors{v}(5:8,k);
        c = vectors{v}(9:12,k);
        d = vectors{v}(13:16,k);
        fourier = inverseFourierDescriptor(a,b,c,d,T,s,ImgSize);        
        red = zeros(ImgSize);
        red(fourier==1) = 255;
        I2 = I;
        I2(:,:,1) = I(:,:,1) + red;
        I2(:,:,2) = I(:,:,2) - red;
        I2(:,:,3) = I(:,:,3) - red;
        figure
        imshow(I2);
        % compute Sørensen-Dice Coefficient between original image and reconstructed 
        dice(v,k) = 2*nnz(canvas&imfill(fourier,'holes'))/(nnz(imfill(fourier,'holes')) + nnz(canvas));     
    end
end
