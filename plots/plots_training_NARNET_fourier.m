%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

% Fourier Descriptors

clc; clear;
cd '/Users/miguelesteras/Desktop/matlab_project/Narnet/'

load('perf_Fourier_Narnet100.mat','record');                                     
Perf = str2num(record);

% Plot performance (validation error during training)
epochs = size(Perf,1);
x = linspace(1,epochs,epochs);
base = Perf(:,2);
y = Perf(:,1);

load('perf_Fourier_Narnet50.mat','record');                                     
Perf = str2num(record);
y2 = Perf(:,1);

figure
plot(x, base, 'Color', [0 0 .5], 'LineWidth', 2); hold on
plot(x, y, 'LineWidth', 3);
plot(x, y2, 'LineWidth', 3);
ylim([0 0.7])
xlim([min(x) max(x)])
grid on
xlabel('Epochs')
ylabel('Root Mean Square Error')
yticks(linspace(0,0.7,8))
set(gca,'fontsize',16)
lg = legend('Baseline','Narnet100','Narnet50');
lg.FontSize = 16;


% Show reconstruction of (numEx) examples  
numEx = 200;
load('/Users/miguelesteras/Desktop/Master Project/data/FourierDescriptorC4.mat','fourierDescriptor');
dataSet = fourierDescriptor;
idx = randi([1 size(dataSet,1)],1,numEx); 
data = dataSet(idx,:);
target = dataSet(idx,end-1);
preFrame = dataSet(idx,end-2);
prediction = [];

% precit and performance
load('net_Fourier_Narnet100.mat','net');
for j = 1:size(data,1)
    query = data(j,1:end-2);
    [xq,xiq,aiq,tq] = preparets(net,{},{},query);
    prediction = [prediction ; net(xq,xiq,aiq)];
end

load('/Users/miguelesteras/Desktop/Master Project/data/cellBodyTraining.mat','cellBodyTraining');
cells = cellBodyTraining(idx,end);
dice = zeros(3,numel(idx));

for k = 1:numel(idx)    
    selection = cells{k};
    selection(:,1) = wrapTo2Pi(selection(:,1) + pi/2); 
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
    
    % inverse fourier transformation
    T = 300;
    s = 0.5;
    ImgSize = size(canvas);
    vectors = {target preFrame prediction}; 
    for v = 1:numel(vectors)        
        a = vectors{v}{k}(1:4);
        b = vectors{v}{k}(5:8);
        c = vectors{v}{k}(9:12);
        d = vectors{v}{k}(13:16);
        BW = inverseFourierDescriptor(a,b,c,d,T,s,ImgSize);
        BW2 = imfill(BW,'holes');
        CC = regionprops(BW2,'BoundingBox');
        CC2 = regionprops(canvas,'BoundingBox');
        BW2area = CC(1).BoundingBox(3)*CC(1).BoundingBox(4);
        canvasArea = CC2(1).BoundingBox(3)*CC2(1).BoundingBox(4);
        ratio = sqrt(canvasArea/BW2area);
        BW3 = imresize(BW2,ratio);
        border = ceil((size(BW3,1)-300)/2);
        BW4 = imcrop(BW3,[border border 299 299]);
%         fourier = bwmorph(BW4,'remove');       
%         color = zeros(ImgSize);
%         color(fourier==1) = 255;
%         I2 = I;
%         if v == 1
%             I2(:,:,1) = I(:,:,1) + color;
%             I2(:,:,2) = I(:,:,2) - color;
%             I2(:,:,3) = I(:,:,3) - color;
%         elseif v ==2
%             I2(:,:,1) = I(:,:,1) - color;
%             I2(:,:,2) = I(:,:,2) - color;
%             I2(:,:,3) = I(:,:,3) + color/2;
%         else 
%             I2(:,:,1) = I(:,:,1) - color;
%             I2(:,:,2) = I(:,:,2) + color/2;
%             I2(:,:,3) = I(:,:,3) - color;
%         end
%         figure
%         imshow(I2);
        % compute Sørensen-Dice Coefficient between original image and reconstructed 
        dice(v,k) = 2*nnz(canvas&BW4)/(nnz(BW4) + nnz(canvas));     
    end
end
save('dice_narnet100_fourier.mat','dice')
