% Stochastic Hidden Markov Models and evaluation of sequence probability
% with forward/backward algorithm
%
% Mikael Mieskolainen, 2010

addpath ./src
close all; clear;


%% Evaluation of sequence probability

% State transition matrix (N x N)
A = [0.9 0.05 0.05;
    0.15 0.8 0.05;
    0.1 0.2 0.7];

% Emission probabilities (each row with sum 1)
B = [0.4 0.3 0.2 0.1;
    0.1 0.7 0.1 0.1;
    0.1 0.2 0.3 0.4];

% Observed sequency
O = [1,1,2,2,2,3,2,1,4,4,4,1,1,1,2,3,2];

% Initial state occupance probabilities
pi = zeros(size(A,1),1); pi(1) = 1;

[P,alfa,beta] = forward_backward(A, B, pi, O);
alfa
beta

