% Stochastic Hidden Markov Models and Viterbi decoding
%
% Mikael Mieskolainen, 2010

addpath ./src
close all; clear;


%% Viterbi decoding of 3 state system

% State transition probabilities
A = [0.500 0.250 0.250;
    0.375 0.125 0.375;
    0.125 0.675 0.375];

% Emission probabilities (each row with sum 1)
B = [0.60 0.20 0.15 0.05;
     0.25 0.25 0.25 0.25;
     0.05 0.10 0.35 0.50];
 
% Observed sequency
O = [1,2,3,4,2,3,4,1,2,3,4,1];

% Initial state occupance probabilities
pi = [0.63, 0.17, 0.20];

% Calculate probability and most likely Viterbi trellis decoded path
[P,q] = viterb(A,B,pi,O)

