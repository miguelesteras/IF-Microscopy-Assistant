%% Phase 5.5 Train Recurrent Network with input delay (pyrenn)
%   ======================================================================
%   Code by Miguel Esteras-Bejar, 07/2017
%   This code is part of the project:
%   'Tracking of temporally occluded or overlapping structures in live cell
%   microscopy'
%   This codes aims to:
%   1. Train a Recurrent network with input delay using the pyrenn
%   package and the fourier descriptors.
%   ======================================================================

rng('default')

inputSiz    = 100;
hidden1Siz  = 500;
hidden2Siz  = 500;
outputSiz   = 100;
dIn         = [0,1,2,3];  % Set of input delays of the neural network
dIntern     = [];  % Set of inernal delays of the neural network
dOut        = [];  % Set of output delays of the neural network

preNet = CreateNN([inputSiz hidden1Siz hidden2Siz outputSiz], dIn, dIntern, dOut);

% Train network with training data, and validate
epochs      = 500;  % max.number of epochs during training to avoid overfitting
terError    = 1e-3; % target prediction error to terminate
net = train_LM(X,Y,preNet,epochs,terError);