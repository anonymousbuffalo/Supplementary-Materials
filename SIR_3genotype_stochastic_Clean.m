function [record, tvec, pA, pB, R0] = SIR_3genotype_stochastic_Clean(mu1, beta, sigma, mu2, r, c, K, outside, sim_time, stoch_pop)
% Gillespie SIR model with 3 host genotypes (AA, AB, BB)
%
% Inputs:
%   mu1      - natural death rate
%   beta     - transmission rate 
%   sigma    - recovery rate 
%   mu2      - disease-induced mortality 
%   r        - baseline fecundity rate 
%   c        - inherent fecundity cost 
%   K        - carrying capacity of the population
%   outside  - external infection pressure
%   sim_time - total simulation time
%   stoch_pop - initial population [S,I,R] for each genotype (3x3)
%
% Outputs:
%   record - matrix recording [time, S1,I1,R1, S2,I2,R2, S3,I3,R3] at each event
%   tvec   - vector of yearly time points
%   pA     - frequency of allele A at each recorded step
%   pB     - frequency of allele B at each recorded step
%   R0     - basic reproduction numbers for each genotype

% Ensuring column vectors
beta    = beta(:);
sigma   = sigma(:);
mu2     = mu2(:);
c       = c(:);
r       = r(:);
outside = outside(:);

% Initialising time and output structures
time = 0;
R0 = beta ./ (sigma + mu1 + mu2); % basic reproduction number
maxsteps = 1.5e7;
record   = zeros(maxsteps, 10);

tvec     = (0:sim_time)'; % yearly time vector
pA       = zeros(length(tvec),1);
pB       = zeros(length(tvec),1);
t_index  = 1;

record(1,:) = [time, reshape(stoch_pop',1,9)];
[pA(t_index), pB(t_index)] = allele_freq(stoch_pop,r,c); 
t_index = t_index + 1;
rec_i = 1;

%% Gillespie Loop
for step = 1:maxsteps
    if time >= sim_time
        break
    end

    S = stoch_pop(:,1); I = stoch_pop(:,2); R = stoch_pop(:,3);
    N = sum(stoch_pop,'all');
    if N == 0
        break
    end

    % Event rates
    lambda       = beta .* (sum(I)/N);
    infection    = lambda .* S;
    outside_inf  = outside .* S;
    recovery     = sigma .* I;
    disease_mort = mu2 .* I;
    deathS       = mu1 .* S;
    deathI       = mu1 .* I;
    deathR       = mu1 .* R;

    % Reproduction
    adult     = r .* (1 - c) .* (S + I + R); 
    B_total   = sum(adult) * max(0, 1 - N/K);

    % Mendelian genotype frequencies
    p = (2*adult(1) + adult(2)) / (2*sum(adult));
    q = 1 - p;
    births = B_total * [p^2; 2*p*q; q^2];

    % Building rate vector (24 events)
    rates = [infection; outside_inf; recovery; disease_mort; ...
             deathS; deathI; deathR; births];

    total_rate = sum(rates);
    if total_rate == 0
        break
    end

    % Gillespie step
    tau = -log(rand) / total_rate;
    time = time + tau;
    r2 = rand()*total_rate;
    ev = find(cumsum(rates) >= r2, 1);

    % Executing events
    if ev <= 3
        g = ev; if stoch_pop(g,1)>0, stoch_pop(g,1)=stoch_pop(g,1)-1; stoch_pop(g,2)=stoch_pop(g,2)+1; end
    elseif ev <= 6
        g = ev-3; if stoch_pop(g,1)>0, stoch_pop(g,1)=stoch_pop(g,1)-1; stoch_pop(g,2)=stoch_pop(g,2)+1; end
    elseif ev <= 9
        g = ev-6; if stoch_pop(g,2)>0, stoch_pop(g,2)=stoch_pop(g,2)-1; stoch_pop(g,3)=stoch_pop(g,3)+1; end
    elseif ev <= 12
        g = ev-9; if stoch_pop(g,2)>0, stoch_pop(g,2)=stoch_pop(g,2)-1; end
    elseif ev <= 15
        g = ev-12; if stoch_pop(g,1)>0, stoch_pop(g,1)=stoch_pop(g,1)-1; end
    elseif ev <= 18
        g = ev-15; if stoch_pop(g,2)>0, stoch_pop(g,2)=stoch_pop(g,2)-1; end
    elseif ev <= 21
        g = ev-18; if stoch_pop(g,3)>0, stoch_pop(g,3)=stoch_pop(g,3)-1; end
    elseif ev == 22
        stoch_pop(1,1) = stoch_pop(1,1)+1;
    elseif ev == 23
        stoch_pop(2,1) = stoch_pop(2,1)+1;
    elseif ev == 24
        stoch_pop(3,1) = stoch_pop(3,1)+1;
    end

    % Recording state
    rec_i = rec_i + 1;
    record(rec_i,:) = [time, reshape(stoch_pop',1,9)];

    % Sampling allele frequencies at yearly intervals
    while t_index <= length(tvec) && time >= tvec(t_index)
        [pA(t_index), pB(t_index)] = allele_freq(stoch_pop,r,c);
        t_index = t_index + 1;
    end
end

record = record(1:rec_i,:);

end
