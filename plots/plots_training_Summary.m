%% Extras. Create Plots
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   ======================================================================

clc; clear;

cd '/Users/miguelesteras/Desktop/Master Project/data'
rhoMatrix = ones(5,200);

% from AE
load('dice_AE.mat','dice')
csvwrite(strcat('dice_AE','.csv'),dice)
rhoMatrix(5,:) = dice(3,:); 

% load dice results for feedforward models
Dicefiles = dir('dice_feedfor*'); 
num_files = length(Dicefiles);
DiceFF = cell(num_files,2);
for i = 1:num_files
    load(Dicefiles(i).name,'dice');                                     
    DiceFF{i,1} = dice;
    DiceFF{i,2} = erase(Dicefiles(i).name,['dice_feedfor_', '.mat']); 
    csvwrite(strcat(DiceFF{i,2},'.csv'),dice)
end

% load dice results for narnet models
Dicefiles = dir('dice_narnet*'); 
num_files = length(Dicefiles);
DiceNarnet = cell(num_files,2);
for i = 1:num_files
    load(Dicefiles(i).name,'dice');                                     
    DiceNarnet{i,1} = dice;
    DiceNarnet{i,2} = erase(Dicefiles(i).name,['dice_narnet', '.mat', '_']); 
    csvwrite(strcat(DiceNarnet{i,2},'.csv'),dice)
end
 save('dice_ff_summary.mat','DiceFF')
 save('dice_narnet.mat','DiceNarnet')
 
