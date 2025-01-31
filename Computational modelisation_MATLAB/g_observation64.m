function [gx] = g_observation64( x, phi, var, inG )
% Currently the following parameters are described as 
% phi is our parameters of interest, our sensitivities = [kR kP kEp kEm kFp kFm];
% var is our variables = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;T-1];
% x is not used in our case as it would be our sensitivities from the last trial
% InG I don't know is function

kR = phi(1);
kP = phi(2);
kEp = phi(3);
kEm = phi(4);
bias = phi(5);


deltaRP = var(1);
deltaE = var(2);
RP_trials = var(3);
Ep_or_Em_trials = var(4);
trial_T = var(5);

% compute equation of our model without evaluation, only observation
% model 1 with interaction param ke kf
SV = (log(1+exp(kR))*(deltaRP)*RP_trials + log(1+exp(kP))*(deltaRP)*(1-RP_trials) - deltaE*Ep_or_Em_trials*(log(1+exp(kEp))) - deltaE*(1-Ep_or_Em_trials)*(log(1+exp(kEm))));
gx = 1/(1+exp(-(SV+bias)));

end

