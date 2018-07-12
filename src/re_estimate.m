% Expectation Maximization Update step for Hidden Markov Models
% 
% The variables as defined in forward_backward_c()
%
% https://en.wikipedia.org/wiki/Baum-Welch_algorithm
%
% Mikael Mieskolainen, SGN-4106/TUT course, 2010

function [A_new, mu_new, sigma2_new] = ...
                                 re_estimate(A, mu, sigma2, O, alfa, beta)

% Re-estimate the transition matrix elements by
A_new = zeros(size(A));
for i = 1:size(A,1)
  for j = 1:size(A,2)
      numm = 0;
      denm = 0;
      for t = 1:length(O)-1
          numm = numm + ...
          alfa(t,i)*A(i,j)*normal_prob(O(t+1),mu(j),sigma2(j))*beta(t+1,j);
          denm = denm + alfa(t,i)*beta(t,i);
      end
      A_new(i,j) = numm/denm;
  end
end

% Re-estimate the means by
mu_new = zeros(size(mu));
for i = 1:length(mu)
    numm = 0;
    denm = 0;
    for t = 1:length(O)
        numm = numm + alfa(t,i)*beta(t,i)*O(t);
        denm = denm + alfa(t,i)*beta(t,i);
    end
    mu_new(i) = numm/denm;
end

% Re-estimate the variances by
sigma2_new = zeros(size(sigma2));
for i = 1:length(sigma2)
    numm = 0;
    denm = 0;
    for t = 1:length(O)
        numm = numm + alfa(t,i)*beta(t,i)*(O(t)-mu_new(i))^2;
        denm = denm + alfa(t,i)*beta(t,i);
    end
    sigma2_new(i) = numm/denm;
end

end