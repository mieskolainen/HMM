% Forward/Backward probability estimation algorithm for Hidden Markov Models
%
% Input:
%   A  = (NxN) Transition probability matrix, A(i,j) ~ prob of transition i->j
%   B  = (NxM) Output emission symbol probabilities,
%        where B(i,j) is the probability of observing the j-th symbol in i-th state
%   pi = Initial state probabilities (Nx1)
%   O  = Observation sequence of size (Tx1) (values between 1 ... M)
%
% Output:
%   P     =  Probability of observation sequence P(O)
%   alfa  =  NxT matrix of alpha values, P(o_1,...,o_t, q_t = i)
%   beta  =  NxT matrix of beta values,  P(o_(t+1),...,o_T | q_t = i)
%
% where
%   T is the observation sequence length
%   N is the number of states
%   M is the number of number of symbols
%
% https://en.wikipedia.org/wiki/Forward-backward_algorithm
%
% Mikael Mieskolainen, SGN-4106/TUT course, 2010

function [P, alfa, beta] = forward_backward(A, B, pi, O)

% Check dimensions
if (size(A,2) ~= size(A,1))
    error('A needs to be a square matrix!');
end
if (size(A,1) ~= size(B,1))
    error('A and B need to have same the same row dimension!');
end
if ((min(O) < 1) || (max(O) > size(B,2)))
    error('Output symbol sequence has values outside the valid range!');
end

% Forward algorithm
[alfa, P] = forward(A, B, pi, O);

% Backward algorithm
[beta] = backward(A, B, O);

end

% ------------------------------------------------------------------------
% Forward algorithm
% ------------------------------------------------------------------------
function [alfa, P] = forward(A, B, pi, O)

    N = size(A,1); % Number of states
    T = length(O); % Observation sequence length
    
    % Initialization
    for i = 1:N
        alfa(1,i) = B(i,O(1))*pi(i);
    end
    
    % Induction
    for t = 1:(T-1)
        for j = 1:N
            summ = 0;
            for i = 1:N
                summ = summ + alfa(t,i)*A(i,j);
            end
            alfa(t+1,j) = summ*B(j,O(t+1));
        end
    end
    
    % Termination
    P = 0;
    for i = 1:N
        P = P + alfa(T,i);        
    end
end


% ------------------------------------------------------------------------
% Backward algorithm
% ------------------------------------------------------------------------
function [beta] = backward(A, B, O)

    N = size(A,1); % Number of state 
    T = length(O); % Observation sequence length
    
    % Initialization
    for i = 1:N
        beta(T,i) = 1;
    end
    
    % Induction
    for t = (T-1):-1:1
        for i = 1:N
            summ = 0;
            for j = 1:N
                summ = summ + A(i,j)*B(j,O(t+1))*beta(t+1,j);
            end
            beta(t,i) = summ;
        end
    end
end
