% Forward/Backward probability estimation algorithm for Hidden Markov
% Models with normal distributed (real valued) outputs
%
% Input:
%   A       =  (N x N) Transition probability matrix,
%              with A(i,j) ~ prob from of transition of state i to j 
%   mu      =  Means of state outputs (N x 1)
%   sigma2  =  Variances of state outputs (N x 1)
%   O       =  Observation sequence (T x 1), real valued
%
% Output:
%   P       =  Probability density value of the sequence O
%   alfa    =  Matrix of alpha values, P(o_1, ... ,o_t, q_t = i)
%   beta    =  Matrix of beta values,  P(o_(t+1), ... ,o_T | q_t = i)
%
% https://en.wikipedia.org/wiki/Forward-backward_algorithm
%
% Mikael Mieskolainen, SGN-4106/TUT course, 2010

function [P, alfa, beta] = forward_backward_c(A, mu, sigma2, pi, O)

% Check dimensions
if (size(A,2) ~= size(A,1))
    error('A needs to be a square matrix!');
end
if ((length(mu) ~= size(A,1) || length(sigma2) ~= size(A,1)))
    error('Mu and sigma2 vectors need to be of size %d x 1 !', N);
end

% Forward algorithm
[alfa, P] = forward(A, mu, sigma2, pi, O);

% Backward algorithm
[beta] = backward(A, mu, sigma2, O);

end

% ------------------------------------------------------------------------
% Forward algorithm
% ------------------------------------------------------------------------
function [alfa, P] = forward(A, mu, sigma2, pi, O)

    N = size(A,1); % Number of states
    T = length(O); % Observation sequence length
    
    % Initialization
    for i = 1:N
        alfa(1,i) = normal_prob(O(1),mu(1),sigma2(1))*pi(i);
    end
    
    % Induction
    for t = 1:(T-1)
        for j = 1:N
            summ = 0;
            for i = 1:N
                summ = summ + alfa(t,i)*A(i,j);
            end
            alfa(t+1,j) = summ*normal_prob(O(t+1),mu(j),sigma2(j));
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
function [beta] = backward(A, mu, sigma2, O)

    N = size(A,1); % Number of states
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
                summ = summ + ...
                    A(i,j)*normal_prob(O(t+1),mu(j),sigma2(j))*beta(t+1,j);
            end
            beta(t,i) = summ;
        end
    end
end
