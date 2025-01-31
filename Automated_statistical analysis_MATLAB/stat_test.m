function [h,p,condition,data1,data2,data3,data1_idx,data2_idx,data3_idx,normal_data] = stat_test(data,split,nb_group,nb_factor,is_paired,allow_transform)
% Code by Arthur Barakat 2022. need publicly available online function(wanova,swtest), personal function (split_data), and toolbox function (boxcox) 
% questions to adress :
% is the data normally distributed ? if no, should I try applying a transformation ?
% How many groups/level (data split in 2 => t-test, 3 or more anova)
% How many factor (features) according to which you split
% Is the data between groups coming from the same subject ? then paired, otherwise unpaired
% does the data has have an identical variance ? otherwise use welch test
% finally, some cases don't have solutions. (if not equal variance, there is no solution). and some
% might not be implemented yet

% Input :
% data, a line vector containing the information you want to test
% split, a struct containing : split.splitter the vector to use to split the data
% split.type, is it a logical vector, or a double ?
% split.threshold, what is the cutoff value on which to split ? if 3 groups, threshold is a vector. 
% nb_group, the number of categories to look into
% nb_factor, the number of features to look into simultaneously
% is_paired, is the data from the two group from the same participant
% allow_transform, in case of none normality, should you transform the data via boxcox ?

% Output :
% h, is the test significant (then h=1)
% p, the p values from the test, in case of anova, p.anova is the value of the anova, and p.post_hoc
% the result of the post hoc test (in case NaN, no post hoc test was performed)
% condition, now everything about the test you did (normal, equal_variance, transformation, etc...)
% data1,data2,data3 are the data split by the split function, if only split in two, then data3 = empty
% data_idx1,2,3 are the index of the split from the original vector containing the data
% normal_data, is equal to data if transformation is not allowed, otherwise is the vector containing
% the transformed version of the data.

h = 0;
p = [];
condition.transform = 'None';
condition.lambda = NaN;
% initialize normal_data as the original dataset. like this a for loop modifies only the correct columns.
normal_data = data;

% split the data in two groups.(data,splitter,condition.type,threshold)
[data1,data2,data3,data1_idx,data2_idx,data3_idx] = split_data(data,split.splitter,split.type,split.threshold);

if 1-isempty(data3)
    data_group = data1_idx+data2_idx*2+data3_idx*3;
else
    data_group = data1_idx+data2_idx*2;
end
is_normal_data1 = 1;
is_normal_data2 = 1;
is_normal_data3 = 1;
condition.is_normal = 1;
condition.is_equal_variance = 0;

if 1-isempty(data1)
    % check if one of the dataset split has a non normal distribution
    is_normal_data1 = 1-swtest(data1);
end
if 1-isempty(data2)
    is_normal_data2 = 1-swtest(data2);
end
if 1-isempty(data3)
    is_normal_data3 = 1-swtest(data3);
end

% if at least one of the dataset has a non normal distribution.
if is_normal_data1 == 0 || is_normal_data2 == 0 || is_normal_data3 == 0
    condition.is_normal = 0;
    if allow_transform == 1
        % since they fail normality, try again by doing a transformation of the original data
        [normal_data,condition.lambda] = boxcox(data');
        % split the data in two groups.(data,splitter,condition.type,threshold)
        [data1,data2,data3,data1_idx,data2_idx,data3_idx] = split_data(normal_data,split.splitter,split.type,split.threshold);
        
        if 1-isempty(data1)
            % check if one of the dataset split has a non normal distribution
            is_normal_data1 = 1-swtest(data1);
        end
        if 1-isempty(data2)
            is_normal_data2 = 1-swtest(data2);
        end
        if 1-isempty(data3)
            is_normal_data3 = 1-swtest(data3);
        end
        % if at least one of the dataset has a non normal distribution, then it fails
        if is_normal_data1 == 0 || is_normal_data2 == 0 || is_normal_data3 == 0
            condition.is_normal = 0;
            % if still fail normality, then use original data, not the transformed one
            [data1,data2,data3,data1_idx,data2_idx,data3_idx] = split_data(normal_data,split.splitter,split.type,split.threshold);
        else
            % transformation worked
            condition.is_normal = 1;
            condition.transform = 'BoxCox';
            condition.lambda = NaN;
        end
    end
end


% now check if the data distribution has equal variance, the equation depends on the nb of groups.
if isempty(data3)
    condition.is_equal_variance = 1-vartest2(data1,data2);
else
    p_val_is_equal_variance = vartestn(data',data_group,'Display','off');
    if p_val_is_equal_variance < 0.05
        condition.is_equal_variance = 1;
    else
        condition.is_equal_variance = 0;
    end
end

% first split if is_normal =1, cannot assume normal distribution
if condition.is_normal == 1
    
    % is there more than 2 groups ?
    if nb_group > 2
        
        % how many factor ?
        if nb_factor >= 2
            % two way anova
            condition.type = 'two_way_anova';
            error('Two way anova not implemented yet')
        else
            % one way anova
            if condition.is_equal_variance == 1
                % condition.type one classical anova, corrected
                condition.type = 'one_way_anova_tuked_correction';
                [p.anova,~,stats_age_blood] = anova1(data,data_group,'off');
                if p.anova <= 0.05
                    h = 1;
                    [tmp,~,~,~] = multcompare(stats_age_blood,'Display','off','CType','tukey-kramer');
                    p.post_hoc = tmp(:,6);
                else
                    p.post_hoc = zeros(nb_group,1);
                end
            else
                % condition.type one classical anova, corrected
                condition.type = 'one_way_welch_anova_no_correction';
                % as wanova doesn't work with Nan. remove the occasional NaN
                data_group(isnan(data)) = [];
                data(isnan(data)) = [];
                [p.anova,~,~,~] = wanova(data',data_group);
                p.post_hoc = NaN(nb_group,1);
            end
        end
        
        % if there is 2 groups
    else
        % is the data paired ?
        if strcmp('paired',is_paired)
            %paired t-test
            condition.type = 'paired_t-test';
            error('Two way anova not implemented yet')
        else
            %unpaired t-test
            if condition.is_equal_variance == 1
                % classical t-test
                condition.type = 'unpaired_t-test';
                [h,p] = ttest2(data1,data2);
            else
                % welch t-test
                condition.type = 'welch_t-test';
                [h,p] = ttest2(data1,data2,'Vartype','unequal');
            end
        end
    end
    % if data not normmaly distributed
else
    % how many groups ?
    if nb_group == 2
        
        % is the data paired or unpaired
        if strcmp('paired',is_paired)
            % wilcoxon signed rank
            condition.type = 'Wilcoxon_signed_rank';
            error('Two way anova not implemented yet')
        else
            % mann whitney U
            condition.type = 'Mann_whitnes_U';
            [p,h] = ranksum(data1,data2);
            
        end
        % if there is 2 groups
    else
        %how many factors ?
        if    nb_factor < 2
            %kruskal wallis
            condition.type = 'Kruskal_wallis_with_tukey_kramer_correction';
            [p.anova,~,stats] = kruskalwallis(data,data_group,'off');
            [tmp,~,~,~] = multcompare(stats,'Display','off','CType','tukey-kramer');
            if p.anova <= 0.05
                h = 1;
                p.post_hoc = tmp(:,6);
            else
                p.post_hoc = zeros(size(tmp,1),1);
            end
        else
            % friedman 2 way
            condition.type = 'Friedman_two_way';
            error('friedman 2 way not implemented yet')
        end
    end
end


end

