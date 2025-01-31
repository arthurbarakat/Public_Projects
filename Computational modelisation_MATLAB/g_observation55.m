function [gx] = g_observation55( x, phi, var, inG )
% Currently the following parameters are described as 
% phi is our parameters of interest, our sensitivities = [kR kP kEp kEm kFp kFm];
% var is our variables = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;T-1;efforts_performed];
% x is used in our case is the sum of efforts performed over time.

% reward sensitivity
kR = phi(1);
%punishment sensitivity
kP = phi(2);
% Physical effort sensitivity
kEp = phi(3);
% Mental effort sensitivity
kEm = phi(4);
%bias for physical task
bias = phi(5);
%physical fatigue
kFp = phi(6);
% mental learning
kLm = phi(7);

% Difference of incentive between the two choices
deltaRP = var(1);
% Difference of effort between the two choices
deltaE = var(2);
% Is it a Reward or Punishment trial as they are intermixed in the 4 blocks
RP_trials = var(3);
% Are we in a physical or mental block
Ep_or_Em_trials = var(4);
% linear increasing value for fatigue
trial_T = var(5);
% we want the sum of effort chosen over time
sum_eff2 = var(9);

% Subjective value computation
SV = (log(1+exp(kR))*(deltaRP)*RP_trials + log(1+exp(kP))*(deltaRP)*(1-RP_trials) - (deltaE)*(log(1+exp(kEp)) + (log(1+exp(kFp))*sum_eff2))*(Ep_or_Em_trials)...
    - (deltaE)*(log(1+exp(kEm)) - (log(1+exp(kLm))*sum_eff2))*(1-Ep_or_Em_trials));
gx = 1/(1+exp(-(SV+bias)));

end
