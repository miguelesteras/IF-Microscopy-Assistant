%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

%% Rho descriptors NARNET performance

load('/Users/miguelesteras/Desktop/Master Project/data/RhoDescriptors.mat','RhoDescriptors');
data = RhoDescriptors;
samples = size(data,1);
maxSam = floor(samples*0.85);

% increasing dataset size
load('/Users/miguelesteras/Desktop/Master Project/models/perf_Rho.mat','rho_performance')
rhoPerf = str2num(rho_performance);
x = [20 50 100 200 400 round(linspace(600,maxSam,11))];
y = rhoPerf;

figure 
plot(x, y(:,1), 'Color', [0.5 0 0], 'LineWidth', 4); hold on
plot(x, y(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x) max(x)])
ylim([0 max(y(:,1))+5])
grid on
xticks([20 50 100 200 400 600 1200 maxSam])
xlabel('Training Data Size')
ylabel('Root Mean Square Error')
set(gca, 'XScale', 'log','fontsize',16)
lg = legend('Model Predictions','Previous Frame');
lg.FontSize = 16;


% increasing epochs
load('/Users/miguelesteras/Desktop/Master Project/models/perf_Rho_1to50Epochs.mat','rho_performance')
rhoPerf = str2num(rho_performance);
x1 = linspace(1,50,50);
y1 = rhoPerf;

figure 
plot(x1, y1(:,1), 'Color', [0.5 0 0], 'LineWidth', 4); hold on
plot(x1, y1(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x1) max(x1)])
ylim([0 max(y1(:,1))+5])
grid on
xticks(linspace(0,50,11));
xlabel('Epochs')
ylabel('Root Mean Square Error')
set(gca,'fontsize',16)
lg = legend('Model Predictions','Previous Frame');
lg.FontSize = 16;


%% Rho descriptors NARNET performance + synthetic data

load('/Users/miguelesteras/Desktop/Master Project/data/RhoDescriptors.mat','RhoDescriptors');
data = RhoDescriptors;
samples = size(data,1);
maxSam = floor(samples*0.85);

% increasing dataset size
load('/Users/miguelesteras/Desktop/Master Project/models/perf_Rho_synthetic.mat','rho_performance')
rhoPerf = str2num(rho_performance);
x = [20 50 100 200 400 round(linspace(600,maxSam,11))];
y2 = rhoPerf;

figure 
plot(x, y2(:,1), 'Color', [.5 .5 0], 'LineWidth', 4); hold on
plot(x, y2(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x) max(x)])
ylim([0 max(y2(:,1))+5])
grid on
xticks([20 50 100 200 400 600 1200 maxSam])
xlabel('Training Data Size')
ylabel('Root Mean Square Error')
set(gca, 'XScale', 'log','fontsize',16)
lg = legend('Model Predictions','Previous Frame');
lg.FontSize = 16;

% increasing epochs
load('/Users/miguelesteras/Desktop/Master Project/models/perf_Rho_1to50Epochs_synthetic.mat','rho_performance')
rhoPerf = str2num(rho_performance);
x1 = linspace(1,50,50);
y3 = rhoPerf;

figure 
plot(x1, y3(:,1), 'Color', [0.5 0.5 0], 'LineWidth', 4); hold on
plot(x1, y3(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x1) max(x1)])
ylim([0 max(y3(:,1))+5])
grid on
xticks(linspace(0,50,11));
xlabel('Epochs')
ylabel('Root Mean Square Error')
set(gca,'fontsize',16)
lg = legend('Model Predictions','Previous Frame');
lg.FontSize = 16;

%% Rho descriptors NARNET performance Combined

figure 
plot(x, y(:,1), 'Color', [0.5 0 0], 'LineWidth', 4); hold on
plot(x, y2(:,1), 'Color', [0.5 0.5 0], 'LineWidth', 4);
plot(x, y(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x) max(x)])
ylim([0 max(y2(:,1))+5])
grid on
xticks([20 50 100 200 400 600 1200 maxSam])
xlabel('Training Data Size')
ylabel('Root Mean Square Error')
set(gca, 'XScale', 'log','fontsize',16)
lg = legend('Original Data (predictions)','Original + Synthetic Data (predictions)','Previous Frame');
lg.FontSize = 16;

figure 
plot(x1, y1(:,1), 'Color', [0.5 0 0], 'LineWidth', 4); hold on
plot(x1, y3(:,1), 'Color', [0.5 0.5 0], 'LineWidth', 4);
plot(x1, y1(:,2), 'Color', [0 0 .5], 'LineWidth', 2);
xlim([min(x1) max(x1)])
ylim([0 max(y3(:,1))+5])
grid on
xticks(linspace(0,50,11));
xlabel('Epochs')
ylabel('Root Mean Square Error')
set(gca,'fontsize',16)
lg = legend('Original Data (predictions)','Original + Synthetic Data (predictions)','Previous Frame');
lg.FontSize = 16;


%% Fourier descriptors NARNET performance

load('/Users/miguelesteras/Desktop/Master Project/models/perf_FourierC4.mat','fourier_performance')
FourierPerf_0_1 = str2num(fourier_performance);

load('/Users/miguelesteras/Desktop/Master Project/models/perf_FourierC4_epochs.mat','fourier_performance')
FourierPerf_1_10 = str2num(fourier_performance);

load('/Users/miguelesteras/Desktop/Master Project/models/perf_FourierC8_epochs.mat','fourier_performance')
FourierPerfC8 = str2num(fourier_performance);

load('/Users/miguelesteras/Desktop/Master Project/models/perf_FourierC10_epochs.mat','fourier_performance')
FourierPerfC10 = str2num(fourier_performance);






%% Boundary descriptors NARNET performance

load('/Users/miguelesteras/Desktop/Master Project/models/boundaryPerf.mat','boundaryPerf')
load('/Users/miguelesteras/Desktop/Master Project/models/net_Boundary_narnet.mat','net')
load('/Users/miguelesteras/Desktop/Master Project/data/BoundaryDescriptors.mat','BoundaryDescriptors');
data = BoundaryDescriptors;
seqLen = size(data,2);
samples = size(data,1);
maxSam = floor(samples*0.85);

x = [50 100 200 400 linspace(600,maxSam,10)]'; % linspace(2*maxSam,10*maxSam,9)];
y = boundaryPerf(3:16,:);
figure 
plot(x,y);
x([6 12],:) = [];
y([6 12],:) = [];

figure 
plot(x, y(:,1), 'Color', [0 0.6 0], 'LineWidth', 4); hold on
plot(x, y(:,2), 'Color', [1 0.5 0.31], 'LineWidth', 2);
xlim([min(x) max(x)])
ylim([5 max(y(:,1))+5])
grid on
xticks([50 100 200 400 600 1000 1600 maxSam])
xlabel('Training Data Size')
ylabel('Root Mean Square Error')
set(gca, 'XScale', 'log','fontsize',16)
lg = legend('Model Predictions','Previous Frame');
lg.FontSize = 16;
