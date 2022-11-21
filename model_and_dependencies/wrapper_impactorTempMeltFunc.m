%% demo for manannan crater simulation in figure 3
clear all; close all;
%% model parameters
% manannan crater key from Cox and Bauer, 2018 simulations
cox_and_bauer_impact_key = '03321';
basal_viscosity = 1e14; % basal viscosity of ice
viscosity_activation_energy = 50e3; % viscosity activation energy
%% model run simulation
impactorTempMeltFunc(cox_and_bauer_impact_key,basal_viscosity, ...
    viscosity_activation_energy)
