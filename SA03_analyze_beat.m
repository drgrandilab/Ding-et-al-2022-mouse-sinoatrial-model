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

[N_trials N_par]=size(all_parameters);

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

if Na_clamp == 1
    disp('Na clamped')
end

duration = 5e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1);

%% Run cycle
% newoutputs = [rr_bpm dVm_max -dVm_min -Vm_min AP_amp...
%         -THR APD APD90 APD50 CL...
%         DD eDD DDR eDDR -MMR...
%         Ca_min Ca_amp Ca_t50 Ca_tau Na_min];
output_names = {'HR', 'UV', '|RR|', '|MDP|', 'AP amp',...
    '|THR|', 'APD', 'APD90', 'APD50', 'CL',...
    'DD', 'EDD', 'DDR', 'lateDDR', '|MRR|',...
    'diast [Ca]', 'CaT amp', 'CaT t50', 'CaT tau', 'diast [Na]'};

output_units = {'bpm', 'mV/ms', 'mV/ms', 'mV', 'mV',...
    'mV', 'ms', 'ms', 'ms', 'ms',...
    'ms', 'ms', 'mV/ms', 'mV/ms', 'mV/ms',...
    'mM', 'mM', 'ms', 'ms', 'mM'};

N_outputs = length(output_names); % number of outputs of beat analysis

all_outputs = zeros(N_trials,N_outputs);

tic
parfor ii=1:N_trials
%for ii=1:100, % plot figure
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    y0n = all_ICs(ii,:);
    par_SA = all_parameters(ii,:); % 19 parameters
    p = [model_index Na_clamp ISO_CCh_flag V_prot input par_block par_SA];
    [t,y] = ode15s(@mouse_SAM_eccODEfile,tspan,y0n,options,p);
    
    currents = calcCurrents(t,y,p,'dVm');
        
    time = t; % (ms)
    Vm = y(:,37); % (mV)
    Ca = y(:,32); % (mM)
    Na = y(:,35); % (mM)
    dVm = currents(:,1); % (mV/ms)
    %figure,plot(t,Vm)

    newoutputs = function_SAN_AP_analysis_single_beat(time,Vm,Ca,Na,dVm,0,2);
    all_outputs(ii,:) = newoutputs;
end

all_outputs
% columns: N outputs
% rows: N trials
toc

%% Saving
%save outputs_matrix_5000_120s_control all_outputs output_names output_units % Control
