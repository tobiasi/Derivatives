%% Input
S     = 100;
K     = 50;
T     = 10;
r     = 0.05;
D     = 0.0;
sigma = 0.15;
type  = 'Call';

% Make the plots sexy
set(groot, 'defaultTextInterpreter', 'LaTex');  
set(groot, 'defaultAxesTickLabelInterpreter', 'LaTex');  
set(groot, 'defaultLegendInterpreter', 'LaTex');  
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);
set(groot, 'defaultAxesLineStyleOrder', '-|:|--');

%% Plot options
optionType = lower(type);
switch lower(optionType)
    case 'call'
        f  = @(x)max(0,x-K*exp(-r*T));
        pT = @(x)max(0,x-K);
        pLb = 1:S*1.5;
        for ST = 1:S*1.5 
            [CallEur(ST), ~] = blsprice(ST,K,r,T,sigma);
            pL(ST) = f(ST);
            pU(ST) = pT(ST); 
        end
        figure(1)
        plot(pL)
        hold on
        plot(CallEur,'k--')
        plot(pU)
        plot(pLb)
        hold off
        legend('Lower bound','Call value' , 'Payoff','Upper bound',...
               'Location','Northwest')
        legend('boxoff') 
        xlim([1 ST])
        ylim([0 ST])
        xticks(K)
        xticklabels({['K = ' int2str(K)]})
        title([type ' option, strike = ' int2str(K)])
        
    case 'put'
        f = @(x)max(0,K*exp(-r*T)-x);
        p  = @(x)max(0,K-x);
        for ST = 1:S*1.5 
            [ ~,PutEur(ST)] = blsprice(ST,K,r,T,sigma);
            pL(ST) = f(ST);
            pU(ST) = p(ST); 
        end
        figure(1)
        plot(pU)
        hold on
        plot(pL)
        plot(PutEur,'k--')
        plot(repmat(K*exp(-r*T),length(PutEur),1))
        hold off
        legend( 'Payoff','Lower bound','Put value','Upper bound','Location','North')
        legend('boxoff') 
        xlim([1 ST])
        ylim([0 ST])
        xticks(K)
        xticklabels({['K = ' int2str(K)]})
        title([type ' option, strike = ' int2str(K)])
        
end 
%% Binomial option price converge to Black-Scholes option price
% Since the Black-Scholes price is a limiting case of the Binomial model as
% n -> \infty, we should also see a clear convergence when increasing the
% number of periods
%==================%
K        = 20;      % Strike
sigma    = 0.30;    % Vol
riskfree = 0.0787;  % risk-free interest rate
divYield = 0;       % Dividend yield
s        = 15;      % Current stock price
type     = 'call';  % Type of option
L        = 1000;    % Amount of intermediate periods
T        = 1;       % Amount of periods
%==================%
start     = 10;
optionVal = nan(L-start+1,1);
for ff = 10:L
    n = ff;
    t = T/n;
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
            payoffnodes(jj-1,indTemp(ii)) = exp(-1*riskfree*t)*...
                       (payoffnodes(jj-1,indTemp(ii-1))*p+...
                        payoffnodes(jj,indTemp(ii-1))*(1-p));
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
figure(2)
plot(ind,oPsmooth)
hold on
plot(blsprices)
hold off
legend('Binomial price', 'Black-Scholes price','Location','North')
legend('boxoff') 
xlim([0 ,L-start])
title(['Binomial option pricing vs. BS - Type = ' , type])
xlabel('Number of periods')

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
figure(4)
plot(dW)
title('Brownian Motion')

% We can also create 1000 of these:
dWmat = cumsum(normrnd(0,sigma_w,nStep,1000));
figure(5)
plot(dWmat)
title('Brownian Motion - 1000 simulations')

% Distribution of final outcomes
figure(6)
histogram(dWmat(end,:))
title('Distribution of final outcomes - BM')

% We see from the histogram that the outcomes appears to be nicely normally
% distributed with a mean of 0. This implies that approximately 50% of the 
% outcomes falls below zero, and is not a nice way to model stock prices 
% which are bounded by zero. A solution is to model the stock prices as a 
% Geometric Brownian Motion (as assumed in the BS-model). The GBM is on the
% form S_t = S_0*exp(mu*t + sigma*W_t)
GBM = exp(mu*t + sigma*dW);
figure(7)
plot(GBM)
title('Geometric Brownian Motion')

% We can also do this 1000 times
GBMmat = exp(mu*t + sigma*dWmat);
figure(8)
plot(GBMmat)
title('Geometric Brownian Motion - 1000 simulations')

% Distribution of final outcomes
figure(9)
histogram(GBMmat(end,:))
title('Distribution of final outcomes - GBM')

% We see from the distribution plot that the outcomes now appears to be
% nicely log-normally distributed, just as we assume the stockprices are


%% 3D plot Greeks
start_t = 0.04;
S1      = 0:0.3:100;
T1      = start_t:0.01:5;
rf      = 0.06;
sigma   = 0.2;
K       = 50;
S       = repmat(S1,length(T1),1);
T       = repmat(T1',1,length(S1));
delta   = nan(length(T1),length(S1)-1);
deltap  = nan(length(T1),length(S1)-1);


for ii = 1:length(T1)
    for jj = 1:length(S1)-1
        [Callprice1, Putprice1] = blsprice(S1(jj),K,rf,T1(ii),sigma); 
        [Callprice2, Putprice2] = blsprice(S1(jj+1),K,rf,T1(ii),sigma);
        f            =  Callprice2-Callprice1;
        fp           =  Putprice2-Putprice1;
        delta(ii,jj) = f;
        deltap(ii,jj)= fp;
    end
end
delta(isnan(delta))   = 0;
deltap(isnan(deltap)) = 0;

figure(10)
mesh(delta)
zlabel('$\frac{dc}{dS}$')
xlabel('$S$')
ylabel('t')
xticks('')
ylabel('$\partial t$')
title(['Delta of a Call option'])


figure(11)
mesh(deltap)
zlabel('$\frac{dp}{dS}$')
xlabel('$S$')
ylabel('t')
xticks('')
ylabel('$\partial t$')
title(['Delta of a Put option'])


for ii = 1:length(T1)
    for jj = 1:length(S1)-2
       gamma(ii,jj)= delta(ii,jj+1)-delta(ii,jj);
       gammap(ii,jj)= deltap(ii,jj+1)-deltap(ii,jj);
    end
end

figure(12)
mesh(gamma)
zlabel('$\frac{d^2c}{dS^2}$')
xlabel('$S$')
ylabel('t')
xticks('')
ylabel('$\partial t$')
title(['Gamma of a Call option'])

figure(13)
mesh(gammap)
zlabel('$\frac{d^2p}{dS^2}$')
xlabel('$S$')
ylabel('t')
xticks('')
ylabel('$\partial t$')
title(['Gamma of a Put option'])


%% Implied volatility
optionprice = 3;
k           = 15;
s           = 17;
r           = 0.03;
guess       = 0.10;
T           = 1;
e           = 0.5;
while abs(e)>0.0001
    
    [n,~] = blsprice(s, k, r, T, guess);
    e     = n - optionprice;
    if e>0
        guess=guess-0.00001;
    else
        guess=guess+0.00001;
    end
    
end

impliedVol = guess;