% main analyze beat
clear
close all
clc

%% Loading initial conditions
% load matrix all_ICs (columns: N state variables, rows: N trials)
load ICs_matrix_5000_120s_control, model_index = 2;
    
% columns: N state variables
% rows: N trials

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
ISO_CCh_flag = 2; % (0 for control, 1 for ISO, 2 for ACh)

V_prot = 0; % 0 for no stimulation
input = 0; % mV, for voltage-clamp protocol

par_block = ones(1,3); % differential block for NKA/NCX/LTCC

par_SA = ones(1,19); % -, for sensitivity analysis

p = [model_index Na_clamp ISO_CCh_flag V_prot input par_block par_SA];

if Na_clamp == 1
    disp('Na clamped')
end

duration = 120e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1);

%% Run cycle
all_ICs_CCh = 0*all_ICs;

tic
parfor ii=1:N_trials
%for ii=1:100, % plot figure
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    y0n = all_ICs(ii,:);
    par_SA = all_parameters(ii,:); % 19 parameters
    p = [model_index Na_clamp ISO_CCh_flag V_prot input par_block par_SA];
    [t,y] = ode15s(@mouse_SAM_eccODEfile,tspan,y0n,options,p);
    %all_ICs(ii,:) = y(end,:);
    [~, mdp_index] = min(y(end-3000:end,37));
    yfinal = y(mdp_index+(length(y(:,37))-3001),:);
    all_ICs_CCh(ii,:) = yfinal;
end

all_ICs = all_ICs_CCh
% columns: N outputs
% rows: N trials
toc

%% Saving
%save ICs_matrix_5000_120s_CCh all_ICs % Control
