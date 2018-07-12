% Stochastic Hidden Markov Models with backward/forward algorithm
% and Viterbi decoding
%
% Mikael Mieskolainen, 2011

addpath ./src
close all; clear;


%% Baum-Welch / EM Maximum Likelihood Training of the system

% Initial estimate of the state transition matrix
A = [0.9 0.05 0.05;
    0.15 0.8 0.05;
    0.1 0.2 0.7];

% Initial estimate of observable mean per state
mu = [1 2 3];

% Initial estimate of observable variance per state
sigma2 = [1 2 3];


% Observation sequence (real valued)
O = [0, 0.5, 1.0, 1.5, 2.0, 1.5, 1.0];

% Initial state occupance probabilities
pi = rand(size(A,1),1); pi = pi / sum(pi);

% Run backward/forward estimation
[P, alfa, beta] = forward_backward_c(A, mu, sigma2, pi, O);

for i = 1:5
    % Re-estimate
    [A, mu, sigma2] = re_estimate(A, mu, sigma2, O, alfa, beta);

    % Calculate the probability of observation again
    [P, alfa, beta] = forward_backward_c(A, mu, sigma2, pi, O);
end

