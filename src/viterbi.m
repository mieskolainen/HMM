% Viterbi algorithm for Hidden Markov Models
%
% Input:
%  A  = (NxN) transition probability matrix, A(i,j) ~ prob of transition i to j
%  B  = output symbol probabilities with
%       B(i,j) being the probability of observing the j-th symbol in i-th state
%  pi = initial probabilities
%  O  = observation sequence of length T (valid elements with 1 ... size(B,2) )
%
% Output:
%  P  = Probability of observing the sequence O
%  q  = Most probable state sequence given the observation
%
% https://en.wikipedia.org/wiki/Viterbi_algorithm
%
% Mikael Mieskolainen, SGN-4106/TUT course, 2010

function [P,q] = viterbi(A, B, pi, O)

T = length(O); % length of observation sequence
N = size(A,1); % number of states

% Initialization
for i = 1:N
    delta(1,i) = pi(i)*B(i,O(1));
    ksi(1,i) = 0;
end

% Recursion
for t = 1:(T-1)
    for j = 1:N
        prods = zeros(1,N);
        for i = 1:N
            prods(i) = delta(t,i)*A(i,j);
        end
        [delta(t+1,j), ksi(t+1,j)] = max(prods);
        delta(t+1,j) = delta(t+1,j)*B(j,O(t+1));
    end
end

% Termination
q = zeros(1,T);
[P, q(T)] = max(delta(T,:));

% Path Back Trace
for t = (T-1):-1:1
   q(t) = ksi(t+1,q(t+1));
end

end
