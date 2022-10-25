% main analysis

clear
close all
clc

%color = [0 0.450980392156863 0.741176470588235]; % BLUE

% To include 2nd population, set this flag to 1:
flag_2nd_pop = 1;

%% Output selection
% 1) rr_bpm 2) dVm_max 3) -dVm_min 4) -Vm_min 5) AP_amp
% 6) -THR 7) APD 8) APD90 9) APD50 10) CL
% 11) DD 12) EDD 13) DDR 14) eDDR 15) -MRR
% 16) Ca_min 17) Ca_amp 18) Ca_t50 19) Ca_tau 20) Na_min

%output_selection=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % all
%output_selection = [2 3 4 5 8 9 10 11 14 16 17 20]; % CL
%output_selection = [1 2 3 4 5 8 9 11 14 16 17 20]; % HR

%output_selection = [1 17 20  2 3 4  5 8 9  11 14 16]; % HR

output_selection = [1 2 4 5 7 20];
% 1) rr_bpm 2) dVm_max 4) -Vm_min 5) AP_amp 7) APD 17) Ca_amp -or- 20) Na_min

%% Loading parameters
load parameter_matrix_5000_0p26 % all_parameters
if flag_2nd_pop == 1
    all_parameters_1p = all_parameters;

    load parameter_matrix_5000_0p26_v2
    all_parameters_2p = all_parameters;

    all_parameters = [all_parameters_1p; all_parameters_2p];
end

[N_trials N_pars] = size(all_parameters);

discard_par = 1;
% Select Parameters
if discard_par == 1
    par_selection = [(1:4) (6:N_pars)];
    all_parameters = all_parameters(:,par_selection);
    parameter_names = parameter_names(par_selection);
    parameter_names{5} = 'GCaL';
    
    [N_trials N_pars] = size(all_parameters);
end

%% Loading outputs
% Control
color_0 = [0 0 0];

load outputs_matrix_5000_120s_control % all_outputs
if flag_2nd_pop == 1
    all_outputs_1p = all_outputs;

    load outputs_matrix_5000_120s_control_v2
    all_outputs_2p = all_outputs;

    all_outputs = [all_outputs_1p; all_outputs_2p];
end

all_outputs = all_outputs(:,output_selection);
all_outputs_0 = all_outputs;

% ISO
color_1 = [255 101 0]/255; % orange

load outputs_matrix_5000_120s_ISO % all_outputs
if flag_2nd_pop == 1
    all_outputs_1p = all_outputs;

    load outputs_matrix_5000_120s_ISO_v2
    all_outputs_2p = all_outputs;

    all_outputs = [all_outputs_1p; all_outputs_2p];
end

all_outputs = all_outputs(:,output_selection);
all_outputs_1 = all_outputs;

% CCh
color_2 = [0 114 189]/255;

load outputs_matrix_5000_120s_CCh % all_outputs
if flag_2nd_pop == 1
    all_outputs_1p = all_outputs;

    load outputs_matrix_5000_120s_CCh_v2
    all_outputs_2p = all_outputs;

    all_outputs = [all_outputs_1p; all_outputs_2p];
end

all_outputs = all_outputs(:,output_selection);
all_outputs_2 = all_outputs;

output_names = output_names(output_selection);
output_units = output_units(output_selection);

N_outputs = length(output_names);
N_figures = ceil(N_outputs/6);

%% Check basic properties
% Control
% Check HR (output(1)), good if > 0
good_trials_HR_0 = (all_outputs_0(:,1) > 0);
% Check APamp (output(4)), good if > 50
good_trials_APamp_0 = (all_outputs_0(:,4) > 50);
% Combine
good_trials_0 = logical(good_trials_HR_0.*good_trials_APamp_0);
% Good count: number of good trials
good_count_0 = sum(good_trials_0)

% ISO
% Check HR (output(1)), good if > 0
good_trials_HR_1 = (all_outputs_1(:,1) > 0);
% Check APamp (output(4)), good if > 50
good_trials_APamp_1 = (all_outputs_1(:,4) > 50);
% Combine
good_trials_1 = logical(good_trials_HR_1.*good_trials_APamp_1);
% Good count: number of good trials
good_count_1 = sum(good_trials_1)

% CCh
% Check HR (output(1)), good if > 0
good_trials_HR_2 = (all_outputs_2(:,1) > 0);
% Check APamp (output(4)), good if > 50
good_trials_APamp_2 = (all_outputs_2(:,4) > 50);
% Combine
good_trials_2 = logical(good_trials_HR_2.*good_trials_APamp_2);
% Good count: number of good trials
good_count_2 = sum(good_trials_2)

%% Combine
good_trials_all = logical(good_trials_HR_0.*good_trials_APamp_0...
    .*good_trials_HR_1.*good_trials_APamp_1.*good_trials_HR_2.*good_trials_APamp_2);

% Good count: number of good trials
good_count_all = sum(good_trials_all)

% Good_parameters: array with parameters from good trials only
good_parameters = all_parameters(good_trials_all,:);

% Good_outputs: array with parameters from good trials only
good_outputs_0 = all_outputs_0(good_trials_all,:);
good_outputs_1 = all_outputs_1(good_trials_all,:);
good_outputs_2 = all_outputs_2(good_trials_all,:);

%% Figures
sp1 = 2; sp2 = 3;
    max_panels = 6;
    
dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs_0(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_0)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Correlation
index = 1;
figure,set(gcf,'color','w')
for i = 1:N_outputs
    for j = 1:N_outputs
        subplot(N_outputs,N_outputs,index)
        set(gca,'box','off','tickdir','out','fontsize',10)
        index = index+1;
        %if i == j
                        
        %else
            plot(good_outputs_0(:,j),good_outputs_0(:,i),'.','Color',color_0)
            if i == N_outputs
                xlabel(output_names{j})
            end
            if j == 1
                ylabel(output_names{i})
            end
        %end 
    
        % Correlation
        % Rule of Thumb for Interpreting the Size of a Correlation Coefficient
        % Size of Correlation   Interpretation
        % .90 to 1.00    		Very high positive (negative) correlation
        % .70 to .90     		High positive (negative) correlation
        % .50 to .70     		Moderate positive (negative) correlation
        % .30 to .50     		Low positive (negative) correlation
        % .00 to .30     		Negligible correlation
        [R, P] = corrcoef(good_outputs_0(:,j),good_outputs_0(:,i));
        %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
        title(['R^2 = ',num2str(R(1,2)^2,3)]);
    end
end

dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs_1(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_1)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Correlation
index = 1;
figure,set(gcf,'color','w')
for i = 1:N_outputs
    for j = 1:N_outputs
        subplot(N_outputs,N_outputs,index)
        set(gca,'box','off','tickdir','out','fontsize',10)
        index = index+1;
        %if i == j
                        
        %else
            plot(good_outputs_1(:,j),good_outputs_1(:,i),'.','Color',color_1)
            if i == N_outputs
                xlabel(output_names{j})
            end
            if j == 1
                ylabel(output_names{i})
            end
        %end 
    
        % Correlation
        % Rule of Thumb for Interpreting the Size of a Correlation Coefficient
        % Size of Correlation   Interpretation
        % .90 to 1.00    		Very high positive (negative) correlation
        % .70 to .90     		High positive (negative) correlation
        % .50 to .70     		Moderate positive (negative) correlation
        % .30 to .50     		Low positive (negative) correlation
        % .00 to .30     		Negligible correlation
        [R, P] = corrcoef(good_outputs_1(:,j),good_outputs_1(:,i));
        %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
        title(['R^2 = ',num2str(R(1,2)^2,3)]);
    end
end

dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs_2(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color_2)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

% Correlation
index = 1;
figure,set(gcf,'color','w')
for i = 1:N_outputs
    for j = 1:N_outputs
        subplot(N_outputs,N_outputs,index)
        set(gca,'box','off','tickdir','out','fontsize',10)
        index = index+1;
        %if i == j
                        
        %else
            plot(good_outputs_2(:,j),good_outputs_2(:,i),'.','Color',color_2)
            if i == N_outputs
                xlabel(output_names{j})
            end
            if j == 1
                ylabel(output_names{i})
            end
        %end 
    
        % Correlation
        % Rule of Thumb for Interpreting the Size of a Correlation Coefficient
        % Size of Correlation   Interpretation
        % .90 to 1.00    		Very high positive (negative) correlation
        % .70 to .90     		High positive (negative) correlation
        % .50 to .70     		Moderate positive (negative) correlation
        % .30 to .50     		Low positive (negative) correlation
        % .00 to .30     		Negligible correlation
        [R, P] = corrcoef(good_outputs_2(:,j),good_outputs_2(:,i));
        %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
        title(['R^2 = ',num2str(R(1,2)^2,3)]);
    end
end

%% Separation (based on ISO and CCh effect)

delta_ISO = good_outputs_1-good_outputs_0;
delta_CCh = good_outputs_2-good_outputs_0;

figure,set(gcf,'color','w')

subplot(2,2,1),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
out_hist = delta_ISO(:,1);
mean_out_hist_1 = mean(out_hist);
std_out_hist_1 = std(out_hist);
histogram(out_hist,'FaceColor',color_1)
xlabel('delta HR w/ ISO (bpm)')
title(['Mean = ',num2str(mean_out_hist_1,3),'; Std = ',num2str(std_out_hist_1,3)])

subplot(2,2,2),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
out_hist = delta_CCh(:,1);
mean_out_hist_2 = mean(out_hist);
std_out_hist_2 = std(out_hist);
histogram(out_hist,'FaceColor',color_2)
xlabel('delta HR w/ CCh (bpm)')
title(['Mean = ',num2str(mean_out_hist_2,3),'; Std = ',num2str(std_out_hist_2,3)])

subplot(2,1,2),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
plot(delta_ISO(:,1),delta_CCh(:,1),'.','Color',color_0)
xlabel('delta HR w/ ISO (bpm)')
ylabel('delta HR w/ CCh (bpm)')
plot([min(delta_ISO(:,1)) max(delta_ISO(:,1))],mean_out_hist_2*[1 1],'r--')
plot(mean_out_hist_1*[1 1],[min(delta_CCh(:,1)) max(delta_CCh(:,1))],'r--')
% Correlation
[R, P] = corrcoef(delta_ISO(:,1),delta_CCh(:,1));
title(['R^2 = ',num2str(R(1,2)^2,3)]);

%% Subset with reduced HR at baseline

mean_HR_ctrl = mean(good_outputs_0(:,1));
std_HR_ctrl = std(good_outputs_0(:,1));

separation = 0.02;

% Control HR
good_trials_HR_ctrl = (good_outputs_0(:,1) < (1-separation)*mean_HR_ctrl);

% ISO effect on HR
good_trials_delta_HR_ISO = (delta_ISO(:,1) > (1+separation)*mean_out_hist_1);


% CCh effect on HR (increased response to CCh)
%good_trials_delta_HR_CCh = (delta_CCh(:,1) < (1+separation)*mean_out_hist_2);

    % CCh effect on HR (decreased response to CCh)
    good_trials_delta_HR_CCh = (delta_CCh(:,1) > (1-separation)*mean_out_hist_2);

    % CCh effect on HR (any response to CCh)
    %good_trials_delta_HR_CCh = (delta_CCh(:,1) > -500);


% Combine ISO/CCh effect
good_trials_ISO_CCh = logical(good_trials_delta_HR_ISO.*good_trials_delta_HR_CCh);
% Good count: number of good trials
good_count_ISO_CCh = sum(good_trials_ISO_CCh)

% Combine
good_trials_mutation = logical(good_trials_HR_ctrl.*good_trials_delta_HR_ISO.*good_trials_delta_HR_CCh);
% Good count: number of good trials
good_count_mutation = sum(good_trials_mutation)

% Good_outputs: array with parameters from good trials only
good_outputs_0_mutation = good_outputs_0(good_trials_mutation,:);
good_outputs_1_mutation = good_outputs_1(good_trials_mutation,:);
good_outputs_2_mutation = good_outputs_2(good_trials_mutation,:);

delta_ISO_mutation = good_outputs_1_mutation-good_outputs_0_mutation;
delta_CCh_mutation = good_outputs_2_mutation-good_outputs_0_mutation;

% Baseline model
baseline_model_FR_control = 408;
baseline_model_FR_ISO = 492;
baseline_model_FR_CCh = 316;

%% Comparison average properties in 2 subgroups

% Mutation
average_0_mutation = mean(good_outputs_0_mutation);
std_0_mutation = std(good_outputs_0_mutation);
average_ISO_mutation = mean(delta_ISO_mutation);
std_ISO_mutation = std(delta_ISO_mutation);
average_CCh_mutation = mean(delta_CCh_mutation);
std_CCh_mutation = std(delta_CCh_mutation);

FR_baseline_mutation = average_0_mutation(1);
FR_baseline_mutation_std = std_0_mutation(1);
deltaFR_ISO_mutation = average_ISO_mutation(1);
deltaFR_ISO_mutation_std = std_ISO_mutation(1);
deltaFR_CCh_mutation = average_CCh_mutation(1);
deltaFR_CCh_mutation_std = std_CCh_mutation(1);

deltaFR_ISO_mutation_rel_array = 100*delta_ISO_mutation(:,1)./good_outputs_0_mutation(:,1);
deltaFR_ISO_mutation_rel = mean(deltaFR_ISO_mutation_rel_array);
deltaFR_ISO_mutation_rel_std = std(deltaFR_ISO_mutation_rel_array);

deltaFR_CCh_mutation_rel_array = 100*delta_CCh_mutation(:,1)./good_outputs_0_mutation(:,1);
deltaFR_CCh_mutation_rel = mean(deltaFR_CCh_mutation_rel_array);
deltaFR_CCh_mutation_rel_std = std(deltaFR_CCh_mutation_rel_array);

control_selection = 0; % set to 1 to restrict control group
if control_selection == 1
    % interval
    delta = 0.20;
    
    % baseline model
    %ref_1 = baseline_model_FR_control; 
    %ref_2 = baseline_model_FR_ISO-baseline_model_FR_control;
    %ref_3 = baseline_model_FR_CCh-baseline_model_FR_control;
    
    % population
    ref_1 = mean_HR_ctrl; 
    ref_2 = mean_out_hist_1;
    ref_3 = mean_out_hist_2;

    % HR
    gtc_delta_HR = (good_outputs_0(:,1) > (1-delta)*ref_1) .* (good_outputs_0(:,1) < (1+delta)*ref_1);
    % ISO effect
    gtc_delta_HR_ISO = (delta_ISO(:,1) > (1-delta)*ref_2) .* (delta_ISO(:,1) < (1+delta)*ref_2);
    % CCh effect
    gtc_delta_HR_CCh = (delta_CCh(:,1) < (1-delta)*ref_3) .* (delta_CCh(:,1) > (1+delta)*ref_3);

    % mean_HR_ctrl mean_out_hist_1 mean_out_hist_2
    
    gtc_logical = logical(gtc_delta_HR.*gtc_delta_HR_ISO.*gtc_delta_HR_CCh);
    good_trials_control = logical((~good_trials_mutation).*gtc_logical);
else
    % Combine - Control group
    good_trials_control = ~good_trials_mutation;
end

% Good count: number of good trials
good_count_control = sum(good_trials_control)

good_outputs_0_control = good_outputs_0(good_trials_control,:);
good_outputs_1_control = good_outputs_1(good_trials_control,:);
good_outputs_2_control = good_outputs_2(good_trials_control,:);

delta_ISO_control = good_outputs_1_control-good_outputs_0_control;
delta_CCh_control = good_outputs_2_control-good_outputs_0_control;

subplot(2,1,2),hold on,
plot(delta_ISO_control(:,1),delta_CCh_control(:,1),'.','Color',color_2)
plot(delta_ISO_mutation(:,1),delta_CCh_mutation(:,1),'.','Color',color_1)
plot(baseline_model_FR_ISO-baseline_model_FR_control,baseline_model_FR_CCh-baseline_model_FR_control,'*')

figure,set(gcf,'color','w')
subplot(1,2,1),hold on
histogram(good_outputs_0_control(:,1),'FaceColor',color_0)
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('HR (bpm)')
ylabel('# models')
title('Control')
subplot(1,2,2),hold on
histogram(good_outputs_0_mutation(:,1),'FaceColor',color_1)
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('HR (bpm)')
ylabel('# models')
title('Mutation')

average_0_control = mean(good_outputs_0_control);
std_0_control = std(good_outputs_0_control);
average_ISO_control = mean(delta_ISO_control);
std_ISO_control = std(delta_ISO_control);
average_CCh_control = mean(delta_CCh_control);
std_CCh_control = std(delta_CCh_control);

FR_baseline_control = average_0_control(1);
FR_baseline_control_std = std_0_control(1);
deltaFR_ISO_control = average_ISO_control(1);
deltaFR_ISO_control_std = std_ISO_control(1);
deltaFR_CCh_control = average_CCh_control(1);
deltaFR_CCh_control_std = std_CCh_control(1);

deltaFR_ISO_control_rel_array = 100*delta_ISO_control(:,1)./good_outputs_0_control(:,1);
deltaFR_ISO_control_rel = mean(deltaFR_ISO_control_rel_array);
deltaFR_ISO_control_rel_std = std(deltaFR_ISO_control_rel_array);

deltaFR_CCh_control_rel_array = 100*delta_CCh_control(:,1)./good_outputs_0_control(:,1);
deltaFR_CCh_control_rel = mean(deltaFR_CCh_control_rel_array);
deltaFR_CCh_control_rel_std = std(deltaFR_CCh_control_rel_array);

figure,set(gcf,'color','w')

subplot(1,3,1),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([FR_baseline_control FR_baseline_mutation])
errorbar([1 2],[FR_baseline_control FR_baseline_mutation],[FR_baseline_control_std FR_baseline_mutation_std],'k.')
ylabel('Baseline HR (bpm)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

subplot(1,3,2),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([deltaFR_ISO_control deltaFR_ISO_mutation])
errorbar([1 2],[deltaFR_ISO_control deltaFR_ISO_mutation],[deltaFR_ISO_control_std deltaFR_ISO_mutation_std],'k.')
ylabel('delta HR w/ ISO (bpm)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

subplot(1,3,3),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([deltaFR_CCh_control deltaFR_CCh_mutation])
errorbar([1 2],[deltaFR_CCh_control deltaFR_CCh_mutation],[deltaFR_CCh_control_std deltaFR_CCh_mutation_std],'k.')
ylabel('delta HR w/ CCh (bpm)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

figure,set(gcf,'color','w')

subplot(1,3,1),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([FR_baseline_control FR_baseline_mutation])
errorbar([1 2],[FR_baseline_control FR_baseline_mutation],[FR_baseline_control_std FR_baseline_mutation_std],'k.')
ylabel('Baseline HR (bpm)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

subplot(1,3,2),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([deltaFR_ISO_control_rel deltaFR_ISO_mutation_rel])
errorbar([1 2],[deltaFR_ISO_control_rel deltaFR_ISO_mutation_rel],[deltaFR_ISO_control_rel_std deltaFR_ISO_mutation_rel_std],'k.')
ylabel('delta HR w/ ISO (%)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

subplot(1,3,3),set(gca,'box','off','tickdir','out','fontsize',12),hold on,
bar([deltaFR_CCh_control_rel deltaFR_CCh_mutation_rel])
errorbar([1 2],[deltaFR_CCh_control_rel deltaFR_CCh_mutation_rel],[deltaFR_CCh_control_rel_std deltaFR_CCh_mutation_rel_std],'k.')
ylabel('delta HR w/ CCh (%)')
set(gca,'XTick',1:2), set(gca,'XTickLabel',{'Ctrl','Mut'}), xlim([0.25 2.75])

%% T-test
% [H1, P1] = ttest2(good_outputs_0_mutation(:,1),good_outputs_0_control(:,1))
% 
% [H2, P2] = ttest2(delta_ISO_mutation(:,1),delta_ISO_control(:,1))
% 
% [H3, P3] = ttest2(delta_CCh_mutation(:,1),delta_CCh_control(:,1))
% 
% [H4, P4] = ttest2(deltaFR_ISO_mutation_rel_array,deltaFR_ISO_control_rel_array)
% 
% [H5, P5] = ttest2(deltaFR_CCh_mutation_rel_array,deltaFR_CCh_control_rel_array)

%% Wilcoxon rank sum test
% p = ranksum(x,y) returns the p-value of a two-sided Wilcoxon rank sum test.
% ranksum tests the null hypothesis that data in x and y are samples from 
% continuous distributions with equal medians, against the alternative that 
% they are not. The test assumes that the two samples are independent. 
% x and y can have different lengths.
% This test is equivalent to a Mann-Whitney U-test.

[P1, H1, stats1] = ranksum(good_outputs_0_mutation(:,1),good_outputs_0_control(:,1));

[P2, H2, stats2] = ranksum(delta_ISO_mutation(:,1),delta_ISO_control(:,1));

[P3, H3, stats3] = ranksum(delta_CCh_mutation(:,1),delta_CCh_control(:,1));

[P4, H4, stats4] = ranksum(deltaFR_ISO_mutation_rel_array,deltaFR_ISO_control_rel_array);

[P5, H5, stats5] = ranksum(deltaFR_CCh_mutation_rel_array,deltaFR_CCh_control_rel_array);

P_array = [P1, P2, P3, P4, P5]
H_array = [H1, H2, H3, H4, H5]

%% Analysis parameteres & boxplot (set flag_boxplot = 1)
% good_parameters parameter_names
% 1) gst 2) gna_ttxs 3) gna_ttxr 4) gcat 5) gcal 
% 6) gh 7) gk1 8) gkr 9) gks 10) gto
% 11) gsus 12) gbna 13) gbca 14) inakmax 15) kNaCa
% 16) ks 17) Pup 18) gkach

group_mutation = find(good_trials_HR_ctrl.*good_trials_delta_HR_ISO.*good_trials_delta_HR_CCh > 0.5);
group_control = find(good_trials_HR_ctrl.*good_trials_delta_HR_ISO.*good_trials_delta_HR_CCh < 0.5);

difference = zeros(1,N_pars);
P_test = zeros(1,N_pars);
H_test = zeros(1,N_pars);

flag_boxplot = 0;
for j = 1:N_pars
	par_perturbations = good_parameters(:,j);
    par_group_mutation = par_perturbations(group_mutation);
    par_group_control = par_perturbations(group_control);
    
    par_group_mutation_mean = mean(par_group_mutation);
    par_group_mutation_std = std(par_group_mutation);
    par_group_control_mean = mean(par_group_control);
    par_group_control_std = std(par_group_control);
    
    if flag_boxplot == 1
        figure,set(gcf,'color','w')
        subplot(1,2,1),boxplot(par_group_control,'Notch','on','Labels','Control','Symbol','ro')
        set(gca,'box','off','tickdir','out','fontsize',12)
        title(['Mean = ',num2str(par_group_control_mean,3),'; Std = ',num2str(par_group_control_std,3)])
        ylim([0.3 3]),ylabel(parameter_names{j})
        subplot(1,2,2),boxplot(par_group_mutation,'Notch','on','Labels','Mutation','Symbol','ro')
        set(gca,'box','off','tickdir','out','fontsize',12)
        title(['Mean = ',num2str(par_group_mutation_mean,3),'; Std = ',num2str(par_group_mutation_std,3)])
        ylim([0.3 3])
    end
    
    difference(j) = par_group_mutation_mean-par_group_control_mean;

    [Pj, Hj] = ranksum(par_group_mutation,par_group_control);
    P_test(j) = Pj; H_test(j) = Hj;
end

figure; set(gcf,'color','w')
subplot(2,1,1)
bar(difference,'FaceColor',color_0)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Difference mean parameteres values (mutation-control)')
set(gca,'XTick',1:N_pars)
set(gca,'XTickLabel',parameter_names)
set(gca,'XLim',[0 N_pars+1])
rotateXLabels( gca(), 90)

subplot(2,1,2)
bar(H_test,'FaceColor',color_0)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Statistical difference? (1 = yes)')
set(gca,'XTick',1:N_pars)
set(gca,'XTickLabel',parameter_names)
set(gca,'XLim',[0 N_pars+1])
rotateXLabels( gca(), 90)
