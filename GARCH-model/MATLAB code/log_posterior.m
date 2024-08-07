function log_posterior = log_posterior(params, data)
    % Extract the parameters
    omega = params(1);
    alpha = params(2);
    beta = params(3);
    if omega <= 0 || alpha <= 0 || beta <= 0 || (alpha + beta) >= 1
        log_posterior = -Inf;
        return;  % No further calculation is necessary
    end

    % GARCH(1,1) log-likelihood calculation as you've defined
    n = length(data);
    sigma2 = zeros(n, 1);
    sigma2(1) = var(data); 
    for t = 2:n
        sigma2(t) = omega + alpha * data(t-1)^2 + beta * sigma2(t-1);
    end
    
    loglikelihood = 0.5*sum(-log(sigma2) - (data.^2) ./ sigma2);
    log_prior_alpha = 0.5*log(alpha)+9*log(1-alpha);
    log_prior_beta = 9*log(beta)+0.5*log(1-beta);
    log_prior = log_prior_alpha + log_prior_beta;

    % Sum of log-likelihood and log-prior
    log_posterior = loglikelihood + log_prior;
end

