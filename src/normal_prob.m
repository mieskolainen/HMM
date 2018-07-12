% Gaussian (normal) density

function P = normal_prob(x, mu, sigma2)
    P = (1/sqrt(2*pi*sigma2))*exp(-(x-mu)^2/(2*sigma2));
end