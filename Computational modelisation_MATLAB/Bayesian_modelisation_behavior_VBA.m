% script to compute behavioral sensitivities and select the best one

%% clean workspace before starting
clc;
clearvars;
close all;
set(0,'defaultfigurecolor',[1 1 1])
%% prepare paths and working directory
% main_folder                 = [pwd]; % you have to be sure that you are in the correct path when you launch the script
results_folder           = % REMOVED FOR PRIVACY
code_folder              = % REMOVED FOR PRIVACY

addpath(results_folder);
addpath(code_folder);
cd(results_folder)

%% set run parameters
% groupe f functions to test
f_fname = {[],[],[],[],[]};
       
% % apply multisession to each block. 
is_multisession = true;

% % groupe g functions to test current best @g_observation55
% test model : with kI, with kI +kE, with kI +kE+bias,with kI,kE+bias+kFlinear, with kI,kE+bias+kFadapted, 
g_fname = {@g_observation66,@g_observation1,@g_observation64,@g_observation65,@g_observation55};

% % Set the number of G parameters. For F param, they are the same number as hidden states in our case
n_G_prm = [2,4,5,7,7];
n_F_prm = [0,0,0,0,0];
n_hiddenStates = [0,0,0,0,0];
n_trialsPerSession = 54;

% set other parameters
nb_functions = length(g_fname);
% change the input date [-2,-1,1,2] to [0,0.25,0.75,1]
rescale_choice = true;
% do you want binary response from the model
binary_answers = false;
% make response binary from model prediction
do_smooth = false;
% should you save the result, and plot it
do_save = true;
do_plot = true;
%% extract behavioral data and put it in a struct
% remove saturators, removes any participants which modeling saturates on 1 run
remove_choice_saturators = true;
remove_choice_saturators_one_run = false;
remove_pred_saturators = false;
% Extract all information for modelling
[p_or_m,nb_participants,all_data,all_files,CID_nb] = extract_behavioral_data(remove_choice_saturators,remove_choice_saturators_one_run,remove_pred_saturators);

%% for each subject and each model, compute modelisation and compute error
RT_mental=NaN(54,72,2);
RT_physical=NaN(54,72,2);
for i_sub = 1:length(all_files)
    i_sub % plot at what subject we currently are.
    [var,all_choices,all_choices_rescaled,deltaRP_idx,deltaRP,deltaE,tmp1,tmp2,all_MVC(i_sub,:),all_NMP(i_sub,:),physical_IP(i_sub)...
        ,mental_IP(i_sub),perf(i_sub),incentive_idx] = reshape_data(i_sub,all_data,rescale_choice,binary_answers,p_or_m);

       % data were acquired in different orders to account for a bias in the experiment, re_order it
    if size(tmp1,1) == 2
    RT_physical(:,i_sub,:) = tmp1';
    elseif size(tmp1,1) == 1
           RT_physical(:,i_sub,1) = tmp1'; 
    end
        if size(tmp2,1) == 2
    RT_mental(:,i_sub,:) = tmp2';
    elseif size(tmp2,1) == 1
           RT_mental(:,i_sub,1) = tmp2; 
    end

  

    % compute decision making reaction time and save it
    RTP(i_sub) = mean(mean(RT_physical(:,i_sub,:)));
    RTM(i_sub) = mean(mean(RT_mental(:,i_sub,:)));

    % save participants choices.
    sub_i_choices(:,i_sub) = all_choices;
    var_p_or_m(:,i_sub) = var(4,:);
    var_saved(:,:,i_sub) = var;
	% choose how many models to test
    for i_model=1:nb_functions
        
        if i_model == 6
                if var_p_or_m(1,i_sub) == 1
                    all_choices(1:54) = NaN;
                    all_choices(109:162) = NaN;
                else
                    all_choices(55:108) = NaN;
                    all_choices(163:216) = NaN;
                end
        end
        % define the options for the VBM toolbox
        [options] = select_options(i_model,all_choices,n_hiddenStates,n_trialsPerSession,n_F_prm,n_G_prm,i_sub,is_multisession);
        
        % since NaN makes the toolbox crash, and that we ignore NaN choices using isYout in the
        % options, we redefine the no choice with a -1 idx, to take it into account later
        all_choices(isnan(all_choices)) = -1;
        var(9,isnan(var(9,:))) = 0;
        var(10,isnan(var(10,:))) = 0;
        var(13,isnan(var(10,:))) = 0;

        % launch toolbox to do bayesian modelisation
        [posterior,out] = VBA_NLStateSpaceModel(all_choices, var, f_fname{i_model}, g_fname{i_model}, options.dim, options);
        

        %put back NaN for easier processing afterwards
        all_choices(all_choices == -1) = NaN;
        
        
        % extract parameters
        sensitivitiesPhi{i_model,i_sub} = posterior.muPhi;
        sensitivitiesTheta{i_model,i_sub} = posterior.muTheta;
        
        % extract free energy matrix and AIC,BIC for model comparison
        free_energy_matrix(i_model, i_sub) = out.F;
        AIC(i_model, i_sub) = out.fit.AIC;
        BIC(i_model, i_sub) = out.fit.BIC;
        R2(i_model, i_sub) = out.fit.R2;
        
        % Extract model predictions for each choice.
        y_hat(:,i_model,i_sub) = out.suffStat.gx;
        
        for i_choice = 1:n_trialsPerSession*4
            % extract predictions as only 2 options.
            if binary_answers == true
                if y_hat(i_choice,i_model,i_sub) >=0 && y_hat(i_choice,i_model,i_sub) < 0.25
                    y_hat(i_choice,i_model,i_sub) = 0;
                elseif y_hat(i_choice,i_model,i_sub) >=0.25 && y_hat(i_choice,i_model,i_sub) < 0.5
                    y_hat(i_choice,i_model,i_sub) = 0;
                elseif y_hat(i_choice,i_model,i_sub) >=0.5 && y_hat(i_choice,i_model,i_sub) <= 0.75
                    y_hat(i_choice,i_model,i_sub) = 1;
                elseif y_hat(i_choice,i_model,i_sub) >0.75 && y_hat(i_choice,i_model,i_sub) <= 1
                    y_hat(i_choice,i_model,i_sub) = 1;
                end
                
            else
                if do_smooth == true
                    % smooth prediction choices in 4 category
                    if y_hat(i_choice,i_model,i_sub) >=0 && y_hat(i_choice,i_model,i_sub) < 0.25
                        y_hat(i_choice,i_model,i_sub) = 0;
                    elseif y_hat(i_choice,i_model,i_sub) >=0.25 && y_hat(i_choice,i_model,i_sub) < 0.5
                        y_hat(i_choice,i_model,i_sub) = 0.25;
                    elseif y_hat(i_choice,i_model,i_sub) >=0.5 && y_hat(i_choice,i_model,i_sub) < 0.75
                        y_hat(i_choice,i_model,i_sub) = 0.75;
                    elseif y_hat(i_choice,i_model,i_sub) >=0.75 && y_hat(i_choice,i_model,i_sub) <= 1
                        y_hat(i_choice,i_model,i_sub) = 1;
                    end
                end
            end
        end
        % compute RMSE or MAE depending on the treatment applied to the data.
        if binary_answers == true
            % warning, make sure all_choices is the rescaled version
            RMSE(i_model,i_sub) = sqrt(nanmean(((all_choices-y_hat(:,i_model,i_sub)').^2)));
            MAE(i_model,i_sub) = nanmean(abs(all_choices-y_hat(:,i_model,i_sub)'));
        else
            
            RMSE(i_model,i_sub) = sqrt(nanmean(((all_choices(~isnan(all_choices))-y_hat(~isnan(all_choices),i_model,i_sub)').^2)));
            MAE(i_model,i_sub) = nanmean(abs(all_choices(~isnan(all_choices))-y_hat(~isnan(all_choices),i_model,i_sub)'));
            
            % Compute physical and mental MAE, perhaps it underfit one system.
            switch p_or_m(i_sub)
                % compute both physical and mental MAE, see where is has a hard time fitting. divide by 12
                % because normally its 216, and it's already grouped (divided by 2) and binned (divided by 9) = 12
                case "p"
                    MAE_physical(i_model,i_sub) = nansum(abs(all_choices(1:54)-y_hat(1:54,i_model,i_sub)')+abs(all_choices(109:162)-y_hat(109:162,i_model,i_sub)'))/108;
                    MAE_mental(i_model,i_sub) = nansum(abs(all_choices(55:108)-y_hat(55:108,i_model,i_sub)')+abs(all_choices(163:216)-y_hat(163:216,i_model,i_sub)'))/108;
                case "m"
                    MAE_mental(i_model,i_sub) = nansum(abs(all_choices(1:54)-y_hat(1:54,i_model,i_sub)')+abs(all_choices(109:162)-y_hat(109:162,i_model,i_sub)'))/108;
                    MAE_physical(i_model,i_sub) = nansum(abs(all_choices(55:108)-y_hat(55:108,i_model,i_sub)')+abs(all_choices(163:216)-y_hat(163:216,i_model,i_sub)'))/108;
            end
        end

        % extract all runs, and run 1 and 2 averaged together
        [nb_ND_physical(i_sub),nb_ND_mental(i_sub),bins_per_trial_choice_physical(:,i_sub),bins_per_trial_choice_mental(:,i_sub),...
        bins_per_trial_pred_physical(:,i_model,i_sub),bins_per_trial_pred_mental(:,i_model,i_sub),grouped_choice_physical_run(:,i_sub),...
        grouped_choice_mental_run(:,i_sub),grouped_pred_physical_run(:,i_model,i_sub),grouped_pred_mental_run(:,i_model,i_sub),...
        saturates_choice_once_per_sub,saturates_choice_twice_per_sub,saturates_pred_per_sub] = extract_individual_sessions(all_choices,sub_i_choices(:,i_sub),y_hat(:,i_model,i_sub),p_or_m(i_sub));
    

        % as parameters are not reset in the struct option, clear it and put new data within
        clear options
        
        % compute participants choices and prediction from the model in heatmaps
        prediction_or_choice = true;
        [phys_choice_mat,ment_choice_mat,nb_trials_choice_mat_physical,nb_trials_choice_mat_mental] = create_choice_mat(y_hat(:,i_model,i_sub)',p_or_m(i_sub),deltaRP_idx,deltaRP,deltaE,prediction_or_choice);
        pred_mat_physical(:,:,i_sub,i_model) = phys_choice_mat;
        pred_mat_mental(:,:,i_sub,i_model) = ment_choice_mat;

    end
    

    
    %Extract choice matrix
    prediction_or_choice = false;
    [tmp_phys_choice_mat,tmp_ment_choice_mat,nb_trials_choice_mat_physical,nb_trials_choice_mat_mental] = create_choice_mat(all_choices,p_or_m(i_sub),deltaRP_idx,deltaRP,deltaE,prediction_or_choice);
    choice_mat_physical(:,:,i_sub) = tmp_phys_choice_mat;
    choice_mat_mental(:,:,i_sub) = tmp_ment_choice_mat;
    nb_trials_mat_physical(:,:,i_sub) = nb_trials_choice_mat_physical;
    nb_trials_mat_mental(:,:,i_sub) = nb_trials_choice_mat_mental;
    
    

end

% Compute confidence for each participant
conf = NaN(216,length(all_files));
for i_sub = 1:length(all_files)
    tmp = sub_i_choices(:,i_sub);
    low_conf = find(tmp == 0.25|tmp == 0.75);
    high_conf = find(tmp == 0|tmp == 1);
    find_NA = find(isnan(tmp));
    conf(low_conf,i_sub) = 0; conf(high_conf,i_sub) = 1; conf(find_NA,i_sub) = NaN;
end
%% compare models
[posterior_BMC, out_BMC] = VBA_groupBMC([free_energy_matrix(:,[1:25,27:69])]);

% define the model of interest for sub analysis
model_i =4;

% prepare perf for plots
[m,p,tmp] = compute_perf_score(perf,incentive_idx,var_saved(2,1:54,:),var_saved(3,1:54,:),sub_i_choices,var_p_or_m,model_i,RT_mental,RT_physical);
%% plot results
if do_plot == true

% plot the results between models
plot_model_selection(out_BMC,AIC,MAE,[],[],BIC,R2);
plot_model_selection_paper1(out_BMC,AIC,MAE,[],[],BIC,R2);

% result for the specific model
plot_choice_and_prediction_matrix(choice_mat_physical,choice_mat_mental,pred_mat_physical,pred_mat_mental,model_i);

plot_choice_prediction_effort_incentive(choice_mat_physical,choice_mat_mental,pred_mat_physical,pred_mat_mental,nb_trials_mat_physical,nb_trials_mat_mental,model_i,CID_nb);

plot_per_trial_and_SV_choice_prediction(grouped_choice_physical_run,grouped_choice_mental_run,grouped_pred_physical_run(:,model_i,:),...
    grouped_pred_mental_run(:,model_i,:),bins_per_trial_choice_physical,bins_per_trial_choice_mental,bins_per_trial_pred_physical(:,model_i,:),...
    bins_per_trial_pred_mental(:,model_i,:),all_NMP(:,1))
end

%% Save all computed datasets
% do you want to save the results ?
if do_save == true
save_data(sensitivitiesPhi,sensitivitiesTheta,CID_nb,p,m,RTP,RTM,all_MVC,all_NMP,physical_IP,mental_IP)

end
cd(results_folder)