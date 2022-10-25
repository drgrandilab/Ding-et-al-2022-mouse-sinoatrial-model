% main PLS

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

output_selection = [1 2 4 5 7 17];
% 1) rr_bpm 2) dVm_max 4) -Vm_min 5) AP_amp 7) APD 17) Ca_amp

%% Loading
plot_index = 0;
% 0 for Control,
% 1 for ISO,
% 2 for CCh

% *************************************************************************
load parameter_matrix_5000_0p26 % all_parameters
if flag_2nd_pop == 1
    all_parameters_1p = all_parameters;

    load parameter_matrix_5000_0p26_v2
    all_parameters_2p = all_parameters;

    all_parameters = [all_parameters_1p; all_parameters_2p];
end

if plot_index == 0 % Control
    load outputs_matrix_5000_120s_control % all_outputs
    if flag_2nd_pop == 1
        all_outputs_1p = all_outputs;
    
        load outputs_matrix_5000_120s_control_v2
        all_outputs_2p = all_outputs;
    
        all_outputs = [all_outputs_1p; all_outputs_2p];
    end

    disp('Control')
    color = [0 0 0];
    discard_par = 1;
end

if plot_index == 1 % ISO
    load outputs_matrix_5000_120s_ISO % all_outputs
    if flag_2nd_pop == 1
        all_outputs_1p = all_outputs;

        load outputs_matrix_5000_120s_ISO_v2
        all_outputs_2p = all_outputs;

        all_outputs = [all_outputs_1p; all_outputs_2p];
    end

    disp('ISO')
    color = [255 101 0]/255; % orange
    discard_par = 1;
end

if plot_index == 2 % CCh
    load outputs_matrix_5000_120s_CCh % all_outputs
    if flag_2nd_pop == 1
        all_outputs_1p = all_outputs;

        load outputs_matrix_5000_120s_CCh_v2
        all_outputs_2p = all_outputs;

        all_outputs = [all_outputs_1p; all_outputs_2p];
    end

    disp('CCh')
    color = [0 114 189]/255; % updated model 1000 - Na Free
    discard_par = 1;
end

% *************************************************************************
disp('--------------------------------')

%% Load parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
% and array 'parameter_names'
% 1) gst 2) gna_ttxs 3) gna_ttxr 4) gcat 5) gcal12 6) gcal13 
% 7) gh 8) gk1 9) gkr 10) gks 11) gto 12) gsus
% 13) gbna 14) gbca 15) inakmax 16) kNaCa 17) ks 18) Pup 19) gkach

%all_parameters(:,end) = all_outputs(:,end)/mean(all_outputs(:,end));

[N_trials N_pars] = size(all_parameters);

%% Load outputs
% load matrix all_outputs (columns: N outputs, rows: N trials)
% and array 'output_names' (and 'output_units')

% %% Outputs - deterministic model
% disp('Outputs with default parameters:')
% disp('--------------------------------')
% 
% disp(['HR = ',num2str(newoutputs(1)), ' bpm'])
% disp(['UV = ',num2str(newoutputs(2)), ' mV/ms'])
% disp(['RR = ',num2str(-newoutputs(3)), ' mV/ms'])
% disp(['MDP = ',num2str(-newoutputs(4)), ' mV'])
% disp(['APamp = ',num2str(newoutputs(5)), ' mV'])
% 
% disp(['THR = ',num2str(-newoutputs(6)), ' mV'])
% disp(['APD = ',num2str(newoutputs(7)), ' ms'])
% disp(['APD90 = ',num2str(newoutputs(8)), ' ms'])
% disp(['APD50 = ',num2str(newoutputs(9)), ' ms'])
% disp(['Cycle length = ',num2str(newoutputs(10)), ' ms'])
% 
% disp(['DD = ',num2str(newoutputs(11)), ' ms'])
% disp(['EDD = ',num2str(newoutputs(12)), ' ms'])
% disp(['DDR = ',num2str(newoutputs(13)), ' mV/ms'])
% disp(['lateDDR = ',num2str(newoutputs(14)), ' mV/ms'])
% disp(['MRR = ',num2str(-newoutputs(15)), ' mV/ms'])
% 
% disp(['diast [Ca]i = ',num2str(1e6*newoutputs(16)), ' nM'])
% disp(['CaT amp = ',num2str(1e6*newoutputs(17)), ' nM'])
% disp(['CaT t50 = ',num2str(newoutputs(18)), ' ms'])
% disp(['CaT tau = ',num2str(newoutputs(19)), ' ms'])
% disp(['diast [Na]i = ',num2str(newoutputs(20)), ' mM'])
% disp('--------------------------------')

%% Outputs - population level
all_outputs_mean = mean(all_outputs);
all_outputs_std_dev = std(all_outputs);

disp('Average outputs at the population level (mean+/-std dev):')
disp('--------------------------------')

disp(['HR = ',num2str(all_outputs_mean(1)),' +/- ',num2str(all_outputs_std_dev(1)),' bpm'])
disp(['UV = ',num2str(all_outputs_mean(2)),' +/- ',num2str(all_outputs_std_dev(2)),' mV/ms'])
disp(['RR = ',num2str(-all_outputs_mean(3)),' +/- ',num2str(all_outputs_std_dev(3)),' mV/ms'])
disp(['MDP = ',num2str(-all_outputs_mean(4)),' +/- ',num2str(all_outputs_std_dev(4)),' mV'])
disp(['APamp = ',num2str(all_outputs_mean(5)),' +/- ',num2str(all_outputs_std_dev(5)),' mV'])

disp(['THR = ',num2str(-all_outputs_mean(6)),' +/- ',num2str(all_outputs_std_dev(6)),' mV'])
disp(['APD = ',num2str(all_outputs_mean(7)),' +/- ',num2str(all_outputs_std_dev(7)),' ms'])
disp(['APD90 = ',num2str(all_outputs_mean(8)),' +/- ',num2str(all_outputs_std_dev(8)),' ms'])
disp(['APD50 = ',num2str(all_outputs_mean(9)),' +/- ',num2str(all_outputs_std_dev(9)),' ms'])
disp(['CL = ',num2str(all_outputs_mean(10)),' +/- ',num2str(all_outputs_std_dev(10)),' ms'])

disp(['DD = ',num2str(all_outputs_mean(11)),' +/- ',num2str(all_outputs_std_dev(11)),' ms'])
disp(['EDD = ',num2str(all_outputs_mean(12)),' +/- ',num2str(all_outputs_std_dev(12)),' ms'])
disp(['DDR = ',num2str(all_outputs_mean(13)),' +/- ',num2str(all_outputs_std_dev(13)),' mV/ms'])
disp(['lateDDR = ',num2str(all_outputs_mean(14)),' +/- ',num2str(all_outputs_std_dev(14)),' mV/ms'])
disp(['MRR = ',num2str(-all_outputs_mean(15)),' +/- ',num2str(all_outputs_std_dev(15)),' mV/ms'])

disp(['diast [Ca]i = ',num2str(1e6*all_outputs_mean(16)),' +/- ',num2str(1e6*all_outputs_std_dev(16)),' nM'])
disp(['CaT amp = ',num2str(1e6*all_outputs_mean(17)),' +/- ',num2str(1e6*all_outputs_std_dev(17)),' nM'])
disp(['CaT t50 = ',num2str(all_outputs_mean(18)),' +/- ',num2str(all_outputs_std_dev(18)),' ms'])
disp(['CaT tau = ',num2str(all_outputs_mean(19)),' +/- ',num2str(all_outputs_std_dev(19)),' ms'])
disp(['diast [Na]i = ',num2str(all_outputs_mean(20)),' +/- ',num2str(all_outputs_std_dev(20)),' mM'])
disp('--------------------------------')

APA = all_outputs(:,5); min_APamp = min(APA)
disp('--------------------------------')

% print_mean = [all_outputs_mean(1:15) 1e6*all_outputs_mean(16) 1e6*all_outputs_mean(17) all_outputs_mean(18:20)]';
% print_std_dev = [all_outputs_std_dev(1:15) 1e6*all_outputs_std_dev(16) 1e6*all_outputs_std_dev(17) all_outputs_std_dev(18:20)]';
% print_both = [print_mean print_std_dev]%;

%% Select outputs
%all_outputs=all_outputs(:,1:end-1*);
%output_names=output_names(1:end-1);
%newoutputs=newoutputs(1:end-1);

all_outputs = all_outputs(:,output_selection);
output_names = output_names(output_selection);
output_units = output_units(output_selection);
%newoutputs=newoutputs(output_selection);

N_outputs = length(output_names);

% Select Parameters
if discard_par == 1
    par_selection = [(1:4) (6:N_pars)];
    all_parameters = all_parameters(:,par_selection);
    parameter_names = parameter_names(par_selection);
    parameter_names{5} = 'GCaL';
    
    [N_trials N_pars] = size(all_parameters);
end

N_figures = ceil(N_outputs/6);

% Istograms
color_pre = [0 0 1];

dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:6
        if dex <= N_outputs
            out_hist = all_outputs(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(2,3,subdex)
            % Plot istogram
            histogram(out_hist,'FaceColor',color_pre)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

index = 1;
figure,set(gcf,'color','w')
for i = 1:N_outputs
    for j = 1:N_outputs
        subplot(N_outputs,N_outputs,index)
        set(gca,'box','off','tickdir','out','fontsize',10)
        index = index+1;
        %if i == j
                        
        %else
            plot(all_outputs(:,j),all_outputs(:,i),'.','Color',color_pre)
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
        [R, P] = corrcoef(all_outputs(:,j),all_outputs(:,i));
        %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
        title(['R^2 = ',num2str(R(1,2)^2,3)]);
    end
end

%% Check basic properties, and define X and Y

% Check
% % Check MDP (output(4)), good if between 59-8 and 59+8 mV
% good_trials_MDP = (all_outputs(:,4)>59-8).*(all_outputs(:,4)<59+8);
% good_trials = logical(good_trials_MDP);
% Check HR (output(1)), good if > 0
good_trials_HR = (all_outputs(:,1) > 0);
% Check APamp (output(4)), good if > 50
good_trials_APamp = (all_outputs(:,4) > 50);
% Combine
good_trials = logical(good_trials_HR.*good_trials_APamp);
% Good count: number of good trials
good_count = sum(good_trials);
% Good_parameters: array with parameters from good trials only
good_parameters = all_parameters(good_trials,:);
X = log(good_parameters);
% Good_outputs: array with parameters from good trials only
good_outputs = all_outputs(good_trials,:);
Y = log(good_outputs);

% % No check
% good_count = N_trials;
% X = log(all_parameters);
% Y = log(all_outputs);

%% Histograms
N_figures = ceil(N_outputs/6);

sp1 = 2; sp2 = 3;
    max_panels = 6;
    
dex = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:sp1*sp2
        if dex <= N_outputs
            out_hist = good_outputs(:,dex);
            mean_out_hist = mean(out_hist);
            std_out_hist = std(out_hist);

            subplot(sp1,sp2,subdex),hold on
            % Plot istogram
            histogram(out_hist,'FaceColor',color)
            set(gca,'box','off','tickdir','out','fontsize',10)
            xlabel([output_names{dex},' (',output_units{dex},')'])
            title(['Mean = ',num2str(mean_out_hist,3),'; Std = ',num2str(std_out_hist,3)])
            dex = dex+1;
        end
    end
end

%% Correlation
index = 1;
figure,set(gcf,'color','w')
for i = 1:N_outputs
    for j = 1:N_outputs
        subplot(N_outputs,N_outputs,index)
        set(gca,'box','off','tickdir','out','fontsize',10)
        index = index+1;
        %if i == j
                        
        %else
            plot(good_outputs(:,j),good_outputs(:,i),'.','Color',color)
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
        [R, P] = corrcoef(good_outputs(:,j),good_outputs(:,i));
        %title(['R^2 = ',num2str(R(1,2)^2,3), '; R = ',num2str(R(1,2),3),'; P = ',num2str(P(1,2),3)]);
        title(['R^2 = ',num2str(R(1,2)^2,3)]);
    end
end

%% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y]=...
    PLS_nipals(X,Y,rank(X));

% % PLS - svds algorithm (2010)
% [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,...
%           Yhat,Yjack,R2x,R2y,RESSy,PRESSy,Q2,r2y_random,rv_random,...
%           Yhat4Press,Yhat4Ress]=PLS_jack_svds(X,Y,rank(X));

N_pars, N_outputs, N_trials, good_count

% Calculate agreement of values predicted by regression (Yhat = Bpls*X) 
% with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

%% Plot
dex1 = 1;
for figdex = 1:N_figures
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex = 1:6
        if dex1 <= N_outputs
            subplot(2,3,subdex)
            % Plot data points
            plot(exp(Y(:,dex1)),exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color);
            xlabel(['Actual ', output_names{dex1}])
            ylabel(['Predicted ', output_names{dex1}])
            title(['R^2 = ',num2str(R2each(dex1),4)])
            set(gca,'box','off','tickdir','out','fontsize',10)
            % Plot identity line
            ylim_ind = get(gca,'ylim') ;
            xlim_ind = get(gca,'xlim') ;
            minpoint = min([ylim_ind(1),xlim_ind(1)]);
            maxpoint = max([ylim_ind(2),xlim_ind(2)]);
            hold on
            plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
            dex1 = dex1+1;
        end
    end
end

figure,set(gcf,'color','w')%,'Position',[50,100,1500,750])
bar(R2each,'FaceColor',color)
set(gca,'box','off','tickdir','out','fontsize',12)
set(gca,'XTick',1:N_outputs)
set(gca,'XTickLabel',output_names)
set(gca,'XLim',[0 N_outputs+1])
ylim([0 1])
rotateXLabels( gca(), 90)
title('R^2 values')

dex2 = 1;
for figdex2 = 1:N_figures,
    figure
    set(gcf,'color','w','Position',[50,100,1500,750])
    for subdex2 = 1:6,
        if dex2 <= N_outputs
            subplot(2,3,subdex2)
            bar(Bpls(:,dex2),'FaceColor',color)
            title(output_names(dex2))
            set(gca,'box','off','tickdir','out','fontsize',10)
            set(gca,'XTick',1:N_pars)
            set(gca,'XTickLabel',parameter_names)
            set(gca,'XLim',[0 N_pars+1])
            rotateXLabels( gca(), 90)
            dex2 = dex2 + 1;
        end
    end
end

% N_figures_p = ceil(N_pars/6);
% dex3 = 1;
% for figdex3 = 1:N_figures_p
%     figure
%     set(gcf,'color','w','Position',[50,100,1500,750])
%     for subdex3 = 1:6
%         if dex3 <= N_pars
%             subplot(2,3,subdex3)
%             bar(Bpls(dex3,:),'FaceColor',color)
%             title(parameter_names(dex3))
%             set(gca,'box','off','tickdir','out','fontsize',10)
%             set(gca,'XTick',1:N_outputs)
%             set(gca,'XTickLabel',output_names)
%             set(gca,'XLim',[0 N_outputs+1])
%             rotateXLabels( gca(), 90)
%             dex3 = dex3 + 1;
%             %ylim([-1 1])
%         end
%     end
% end

%% Plot X, B and Y matrices
plot_matrix = 0;

if plot_matrix == 1
    figure; set(gcf,'color','w')
    imagesc(all_parameters); colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Parameters (X)');
    xlabel('Parameters');
    ylabel('Trials');
    %set(gca,'YTick',(1:N_trials))
    set(gca,'XTick',(1:N_pars))
    set(gca,'XTickLabel',parameter_names)
    rotateXLabels( gca(), 90)
    colorbar

    figure; set(gcf,'color','w')
    imagesc(Bpls'); colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Regression coefficients (B)');
    ylabel('Parameters');
    xlabel('Outputs');
    set(gca,'YTick',(1:N_outputs))
    set(gca,'YTickLabel',output_names)
    set(gca,'XTickLabel',parameter_names)
    set(gca,'XTick',(1:N_pars))
    rotateXLabels( gca(), 90)
    colorbar

    norm_outputs=all_outputs;
    for ii=1:N_trials
        norm_outputs(ii,:)=all_outputs(ii,:)./newoutputs;
    end

    figure; set(gcf,'color','w')
    %imagesc(all_outputs); colormap jet;
    imagesc(norm_outputs); colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Outputs (Y)');
    xlabel('Outputs');
    ylabel('Trials');
    %set(gca,'YTick',(1:N_trials))
    set(gca,'XTick',(1:N_outputs))
    set(gca,'XTickLabel',output_names)
    rotateXLabels( gca(), 90)
    colorbar
end

%% Plot perturbation effect
plot_effect = 0;

if plot_effect == 1,
    var_index = 5; % Na
    %var_value = newoutputs(var_index);
    var_value = 1; % for normalized values

    if discard_par == 1,
        par_indexes = [15-1 6-1 9-1 13-1 4];
    else
        par_indexes = [15 6 9 13 4];
    end

    par_index = par_indexes(1); % 15; % NKA
    B = Bpls(par_index,var_index);

    scale = (1-0.30:0.01:1+0.30);
    modulation = var_value.*scale.^B;

    figure; set(gcf,'color','w')
    hold on, plot(scale,modulation)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Scale factor (-)');
    %ylabel('Output variation (-)');
    ylabel('APD (ms)');

    par_index_2 = par_indexes(2); % 6; % ICaL
    B2 = Bpls(par_index_2,var_index)
    mod2 = var_value.*scale.^B2;

    par_index_3 = par_indexes(3); % 9; % IKr
    B3 = Bpls(par_index_3,var_index)
    mod3 = var_value.*scale.^B3;

    par_index_4 = par_indexes(4); % 13; % IbkgNa
    B4 = Bpls(par_index_4,var_index)
    mod4 = var_value.*scale.^B4;
    
    par_index_5 = par_indexes(5); % 4; % ICaT
    B5 = Bpls(par_index_5,var_index)
    mod5 = var_value.*scale.^B5;

    plot(scale,mod2,scale,mod3,scale,mod4,scale,mod5)

    %legend('vNKA','GCaL','GKr','GbkgNa')
    legend(parameter_names(par_indexes))
    xlim([scale(1) scale(end)])
end

%%
plot_effect_2 = 0;

if plot_effect_2 == 1,
    var_index = 1; % FR
    %var_value = newoutputs(var_index);
    var_value = 1; % for normalized values

    if discard_par == 1,
        par_indexes = [15-1 6-1 9-1 13-1 4];
    else
        par_indexes = [15 6 9 13 4];
    end

    par_index = par_indexes(1); % 15; % NKA
    B = Bpls(par_index,var_index);

    scale = (1-0.30:0.01:1+0.30);
    modulation = var_value.*scale.^B;

    figure; set(gcf,'color','w')
    hold on, plot(scale,modulation)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Scale factor (-)');
    %ylabel('Output variation (-)');
    ylabel('FR (bpm)');

    par_index_2 = par_indexes(2); % 6; % ICaL
    B2 = Bpls(par_index_2,var_index)
    mod2 = var_value.*scale.^B2;

    par_index_3 = par_indexes(3); % 9; % IKr
    B3 = Bpls(par_index_3,var_index)
    mod3 = var_value.*scale.^B3;

    par_index_4 = par_indexes(4); % 13; % IbkgNa
    B4 = Bpls(par_index_4,var_index)
    mod4 = var_value.*scale.^B4;
    
    par_index_5 = par_indexes(5); % 4; % ICaT
    B5 = Bpls(par_index_5,var_index)
    mod5 = var_value.*scale.^B5;

    plot(scale,mod2,scale,mod3,scale,mod4,scale,mod5)

    %legend('vNKA','GCaL','GKr','GbkgNa')
    legend(parameter_names(par_indexes))
    xlim([scale(1) scale(end)])
end