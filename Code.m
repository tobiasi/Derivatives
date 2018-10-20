%% Input
S     = 100;
K     = 80;
T     = 10;
r     = 0.05;
D     = 0.0;
sigma = 0.15;
type  = 'Put';

%% Plot options
optionType = lower(type);
switch lower(optionType)
    case 'call'
        f = @(x)max(0,x-K*exp(-r*T));
    case 'put'
        f = @(x)max(0,K*exp(-r*T)-x);
        p = @(x)max(0,K-x);
end 
for ST = 1:S*1.5 % Adjust 1.5 if you want more along the x-axis
    [CallEur(ST),PutEur(ST)] = blsprice(ST,K,r,T,sigma);
    po1(ST) = f(ST);
    po2(ST) = p(ST); 
end
prices    = [CallEur',PutEur'];
c         = {'call','put'};
ind       = find(ismember(c,optionType));
pricesOut = prices(:,ind);


set(groot, 'defaultTextInterpreter', 'LaTex');  
set(groot, 'defaultAxesTickLabelInterpreter', 'LaTex');  
set(groot, 'defaultLegendInterpreter', 'LaTex');  
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);
set(groot, 'defaultAxesLineStyleOrder', '-|:|--');
plot(po1)
hold on
plot(po2)
plot(pricesOut,'k--')
hold off
legend('Lower bound', 'Payoff', [type,' value'],'Location','North')
legend('boxoff') 
xlim([1 ST])
ylim([0 ST])
xticks(K)
xticklabels({['K = ' int2str(K)]})
title([type ' option, strike = ' K])


%% Binomial option price converge to Black-Scholes option price
% Since the Black-Scholes price is a limiting case of the Binomial model as
% n -> \infty, we should also see a clear convergence when increasing the
% number of periods
L         = 1000;
start     = 10;
optionVal = nan(L-start+1,1);
for ff = 10:L
    % --------------------------------------- %
    n        = ff;
    T        = 1;
    t        = T/n;
    K        = 20;
    sigma    = 0.30;
    riskfree = 0.0787;
    divYield = 0;
    s        = 15;
    type     = 'call';
    % ---------------------------------------- %
    switch type
        case 'call'
            f = @(x)max(0,x-K);
        case 'put'
            f = @(x)max(0,K-x);
    end 
    u       = exp((riskfree-divYield)*t + sigma*t^.5);
    d       = exp((riskfree-divYield)*t - sigma*t^.5);
    mat     = nan(101,101);
    nodes   = nan(n+1,n+1);
    indTemp = flip(1:1:n+1);
    nodes   = u.^(n+1:-1:1).*d.^(1:n+1);
    p       = (exp((riskfree-divYield)*t)-d)/(u-d);
    payoffnodes        = nan(n+1,n+1);
    payoffnodes(:,end) = max(0,nodes(:,end)*s-K);
    payoffnodes(:,end) = f(nodes*s);
    
    for ii = 2:length(nodes)
        for jj = 2:indTemp(ii)+1
            payoffnodes(jj-1,indTemp(ii)) = exp(-1*riskfree*t)*(payoffnodes(...
                    jj-1,indTemp(ii-1))*p+payoffnodes(jj,indTemp(ii-1))*(1-p));
        end
    end
    
    optionVal(ff-9) = payoffnodes(1);


end

ind            = 1:10:L;
optVal2plot    = optionVal(ind);
p              = polyfit(transpose(1:100),optVal2plot,7);
oPsmooth       = polyval(p,1:100);
[bsCall,bsPut] = blsprice(s,K,riskfree,T,sigma);
c              = {'call','put'};
indT           = find(ismember(c,type));
prices         = [bsCall,bsPut];
blsprices      = repmat(prices(:,indT),L-start+1,1);
plot(ind,oPsmooth)
hold on
plot(blprice)
hold off
legend('Binomial price', 'Black-Scholes price','Location','North')
legend('boxoff') 
xlim([0 ,L-start])
title('Binomial option pricing vs. Black-Scholes')
xlabel('Number of periods')


%% Brownian Motion (Wiener process)
% Assuming that S_t ~ GBM with unconditional variance
% Set parameters
nStep   = 10;
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













































