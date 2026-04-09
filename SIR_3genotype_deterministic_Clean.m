function [T,Y, pA, pB, R0] = SIR_3genotype_deterministic_Clean(mu1,beta,sigma,mu2, r, c, K, time, det_pop)
% Basic Deterministic SIR model with 3 host genotypes (AA, AB, BB)
% Including the following parameters
% - natural death rate
% - transmission rate
% - recovery rate
% - disease-induced mortality
% - Baseline fecundity rate 
% - Inherent fecundity cost
% - Time
% - Carrying capacity
% This model also simulates Mendelian inheritance of genotypes 

% Packing parameters into structure
pars.mu1 = mu1;          
pars.beta = beta;       
pars.sigma = sigma;     
pars.mu2 = mu2;         
pars.r = r;             
pars.c = c;             
pars.K = K;             

% Computing basic reproduction numbers (R0) for each genotype
R0 = beta ./ (sigma + mu1 + mu2);

% Checking initial population
if numel(det_pop) ~= 9
    error('init_pop must be a 9-element vector: [S_AA I_AA R_AA S_AB I_AB R_AB S_BB I_BB R_BB]');
end

% Integrating system
[T, Y] = ode15s(@equations, 0:1:time, det_pop, ...
    odeset('nonnegative', 1:9, 'RelTol', 1e-6, 'AbsTol', 1e-9), ...
    pars);

% Computing allele frequencies over time
N_AA = Y(:,1) + Y(:,2) + Y(:,3);
N_AB = Y(:,4) + Y(:,5) + Y(:,6);
N_BB = Y(:,7) + Y(:,8) + Y(:,9);
N_tot = N_AA + N_AB + N_BB;

pA = (2*N_AA + N_AB) ./ (2*N_tot);
pB = 1 - pA;

end


function dydt = equations(~, y, pars)
% Unpacking parameters 
mu1 = pars.mu1;      
beta = pars.beta;     
sigma = pars.sigma;    
mu2 = pars.mu2;    
r = pars.r;           
c = pars.c;           
K = pars.K;           

% Current state 
S = [y(1), y(4), y(7)];
I = [y(2), y(5), y(8)];
R = [y(3), y(6), y(9)];
N_i = S + I + R;      
N = sum(N_i);         
if N == 0
    dydt = zeros(9,1);
    return
end

% Force of infection
lambda_i = beta .* (sum(I)/N); 

% Density-dependent birth rates and genotype-specific adult contributions
adult_contrib = r .* (1 - c) .* (S + I + R);  
B_total = sum(adult_contrib) * max(0, (1 - N / K));  

% Mendelian inheritance
p = (2*adult_contrib(1) + adult_contrib(2)) / (2*sum(adult_contrib));
q = 1 - p;
f_offspring = [p^2, 2*p*q, q^2];
births = B_total .* f_offspring; 

% Differential equations
dydt = zeros(9,1);
for i = 1:3
    Si = S(i);
    Ii = I(i);
    Ri = R(i);
    
    dydt(3*(i-1)+1) = births(i) - (lambda_i(i) + mu1)*Si;            % dS/dt
    dydt(3*(i-1)+2) = lambda_i(i)*Si - (sigma(i) + mu1 + mu2(i))*Ii; % dI/dt
    dydt(3*(i-1)+3) = sigma(i)*Ii - mu1*Ri;                           % dR/dt
end
end