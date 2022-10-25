% main generates random parameters
clear
close all
clc

%% Parameters
% 1) gst 2) gna_ttxs 3) gna_ttxr 4) gcat 5) gcal12 6) gcal13 
% 7) gh 8) gk1 9) gkr 10) gks 11) gto 12) gsus
% 13) gbna 14) gbca 15) inakmax 16) kNaCa 17) ks 18) Pup 19) gkach
parameter_names = {'Gst','GNa1.1','GNa1.5','GCaT','GCaL1.2','GCaL1.3',...
    'Gf','GK1','GKr','GKs','Gto','Gsus',...
    'GbNa','GbCa','vNKA','vNCX','vRyR','vSERCA','GKACh'} ;

n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variation
variations = 5000; % number of trials

%sigmaG = 0.1*ones(1,n_parameters); % standard deviation for parameters
sigmaG = 0.26*ones(1,n_parameters); % standard deviation for parameters
%sigmaG = 0.3*ones(1,n_parameters); % standard deviation for parameters

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end

all_parameters %size(all_parameters)
% columns: N parameters
% rows: N trials

%% Saving
%save parameter_matrix_5000_0p26 all_parameters parameter_names
