% main obtain ICs
clear
close all
clc

%% Loading initial conditions
load yfin_Kharche_optimized, model_index = 2;  % model_index = 2;
    
y0n = yfinal;
N_state_vars = length(y0n);

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load parameter_matrix_5000_0p26 % sigma 0.26

[N_trials N_par] = size(all_parameters);

%% Input parameters
Na_clamp = 0; % [0 for free Na, 1 for Na clamp]
if Na_clamp == 1
    disp('Na clamped')
end

% Isoproterenol/Carbachol administration
ISO_CCh_flag = 0; % (0 for control, 1 for ISO, 2 for ACh)

V_prot = 0; % 0 for no stimulation
input = 0; % mV, for voltage-clamp protocol

par_block = ones(1,3); % differential block for NKA/NCX/LTCC

par_SA = ones(1,19); % -, for sensitivity analysis

p = [model_index Na_clamp ISO_CCh_flag V_prot input par_block par_SA];

duration = 120e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1);

%% Run cycle
all_ICs = zeros(N_trials,N_state_vars);

% tic
% for ii=1:N_trials,
%     X = sprintf('Run %d on %d',ii,N_trials); disp(X)
%     par_SA = all_parameters(ii,:); % 18 parameters
%     p = [model_index Na_clamp DB OLD V_prot input par_SA];
%     [t,y] = ode15s(@KharcheSAN_eccODEfile,tspan,y0n,options,p);
%     all_ICs(ii,:) = y(end,:);
%     
%     figure,
%     subplot(4,1,1),plot(t,y(:,37)),ylabel('Em')
%     subplot(4,1,2),plot(t,y(:,32)),ylabel('Ca')
%     subplot(4,1,3),plot(t,y(:,33)),ylabel('Ca rel')
%     subplot(4,1,4),plot(t,y(:,35)),ylabel('Na')
% end

tic
parfor ii=1:N_trials
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    par_SA = all_parameters(ii,:); % 19 parameters
    p = [model_index Na_clamp ISO_CCh_flag V_prot input par_block par_SA];
    [t,y] = ode15s(@mouse_SAM_eccODEfile,tspan,y0n,options,p);
    %all_ICs(ii,:) = y(end,:);
    [~, mdp_index] = min(y(end-3000:end,37));
    yfinal = y(mdp_index+(length(y(:,37))-3001),:);
    all_ICs(ii,:) = yfinal;
end

all_ICs
% columns: N state variables
% rows: N trials
toc

%% Saving
%save ICs_matrix_5000_120s_control all_ICs % Control
