%% Input
S          = 100;
K          = 80;
T          = 10;
r          = 0.05;
D          = 0.0;
sigma      = 0.15;
optionType = 'call';

%% Plot options
    switch optionType
        case 'call'
            f = @(x)max(0,x-K*exp(-r*T));
        case 'put'
            f = @(x)max(0,K*exp(-r*T)-x);  
    end 
for ST = 1:S*1.5 % Adjust 1.5 if you want more along the x-axis
    [CallEur(ST),PutEur(ST)] = blsprice(ST,K,r,T,sigma);
    PayoffEur(ST)            = f(ST);
end
prices    = [CallEur',PutEur'];
c         = {'call','put'};
ind       = find(ismember(c,optionType));
pricesOut = prices(:,ind);

plot(PayoffEur)
hold on
plot(pricesOut,'k--')
hold off
xlim([1 ST])
ylim([0 ST])
xticks([K])
xticklabels({['K = ' int2str(K)]})


%% Brownian Motion (Wiener process)
% Assuming that S_t ~ GBM with unconditional variance
% Set parameters
nStep   = 1000;
T       = 1;
mu      = .15;
sigma_w = sqrt(T/nStep);
t       = (1:nStep)'./nStep;

% A brownian motion is simply a process which is the sum of a stochastic
% process defined on a probability space which fullfills three properties: 
% i)   Integrability, |E[W_T]| < +\infty  
% ii)  Adaptive to filtration, F_t. 
% iii) Martingale, E[W_t|F_s] = W_s for all t>=s. 
%
% We can create this process by simply take the cumulative sum of a
% symmetric stochastic process with mean zero (because of symmetry) and
% variance sigma_w.
dW = cumsum(normrnd(0,sigma_w,nStep,1));
plot(dW)

% We can also create 1000 of these:
dWmat = cumsum(normrnd(0,sigma_w,nStep,1000));
plot(dWmat)

% Distribution of final outcomes
histogram(dWmat(end,:))

% We see from the histogram that the outcomes appears to be nicely normally
% distributed with a mean of 0. This implies that approximately 50% of the 
% outcomes falls below zero, and is not a nice way to model stock prices 
% which are bounded by zero. A solution is to model the stock prices as a 
% Geometric Brownian Motion (as assumed in the BS-model). The GBM is on the
% form S_t = S_0*exp(mu*t + sigma*W_t)
GBM = exp(mu*t + sigma*dW);
plot(GBM)

% We can also do this 1000 times
GBMmat = exp(mu*t + sigma*dWmat);

% Distribution of final outcomes
histogram(GBMmat(end,:))

% We see from the distribution plot that the outcomes now appears to be
% nicely log-normally distributed, just as we assume the stockprices are















































