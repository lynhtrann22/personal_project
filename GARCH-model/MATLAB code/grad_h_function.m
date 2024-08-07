function [grad_h_lambda, h_lambda] = grad_h_function(theta,y,mu,L,n)
% compute the gradient of h_lambda function w.r.t. coefficient vector beta.
% Also returns the value of h_lambda function
d = length(theta);
Sigma = L*(L');

w = exp (theta(1));
alpha = exp(theta(2))*exp(theta(3))/(1 + exp(theta(2)) + exp(theta(3)) + exp(theta(2))*exp(theta(3)));
beta = exp(theta(2))/(1 + exp(theta(2)) + exp(theta(3)) + exp(theta(2))*exp(theta(3)));

sigma2_der = zeros(n,d);
sigma2 = zeros(n,1);
sigma2(1) = y(1);

aux1 = exp(-theta(2))/(1 + exp(-theta(2)));
aux2 = exp(-theta(2))/(1 + exp(-theta(2)));
aux3 = exp(-theta(3))/(1 + exp(-theta(3)));
aux4 = - exp(theta(3))/(1 + exp(theta(3))); 

for t = 2:n
    sigma2(t) = w + alpha.*y(t-1).^2 + beta.*sigma2(t-1);
    sigma2_der(t,1) = w + beta*sigma2_der(t-1,1);
    sigma2_der(t,2) = alpha*aux1*y(t-1).^2 + beta*aux2*sigma2(t-1) + beta*sigma2_der(t-1,2);
    sigma2_der(t,3) = alpha*aux3*y(t-1).^2 + beta*aux4*sigma2(t-1) + beta*sigma2_der(t-1,3);
end

log_prior = 1; 
log_llh = -n/2*log(2*pi) - 1/2*sum(log(sigma2)) -1/2*sum((y.^2)./sigma2);
log_q = -d/2.*log(2*pi) + log(det(L)) - (1/2)*(theta - mu)'*Sigma*(theta - mu);
h_lambda = log_prior + log_llh - log_q;

sigma2_matrix = sigma2.*ones(n,3);
y_matrix = y.*ones(n,3);
grad_log_llh = 1/2*(sum(1./sigma2_matrix.*sigma2_der.*(y_matrix.^2./sigma2_matrix - 1)));
grad_log_q = -Sigma*(theta - mu);
grad_log_prior = 0; 
grad_h_lambda = grad_log_prior+grad_log_llh'-grad_log_q;

end
    