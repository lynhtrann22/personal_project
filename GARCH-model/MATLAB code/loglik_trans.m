function loglik_trans = loglik_trans(params,data)
    % Extract the parameters
    theta1 = params(1);
    theta2 = params(2);
    theta3 = params(3);

    % Transform said params into original params
    omega = exp(theta1);
    alpha = exp(theta2)*exp(theta3)/(1+exp(theta2)+exp(theta3)+exp(theta2)*exp(theta3));
    beta = exp(theta2)/(1+exp(theta2)+exp(theta3)+exp(theta2)*exp(theta3));

    % Compute sigma2_t (variance of the return series)
    n = length(data);
    sigma2 = zeros(n, 1);
    sigma2(1) = var(data); % assuming you want to start with some kind of sensible number
    for t = 2:n
        sigma2(t) = omega + alpha * data(t-1)^2 + beta * sigma2(t-1);
    end

    % Compute loglik
    loglik_trans = 0.5*sum(-log(sigma2) - (data.^2) ./ sigma2);


