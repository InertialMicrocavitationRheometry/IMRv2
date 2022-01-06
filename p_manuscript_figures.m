%Main program to generate and save the figures
% clearing workspace, closing figures, and clear command window
clear; close all; clc;

% figure 1 is the data fitting for the Carreau viscosity parameters using
% data from Biro (1982), Merrill (1969), and Skalak et al. (1981)
% figures are saved in folder: fitfigs
run m_BMS_Carreaufit.m

% figure 2 is the data fitting for the Carreau viscosity parameters using
% the Kim et al. 2017 (Sensors) data for three hematocrit levels
% figures are saved in folder: fitfigs
run m_Kim_Carreaufit.m

% figures are saved in folder: nonnewt_length
run m_nonnewt_length.m

% figures are saved in folder: shear_analysis
run m_shear_analysis.m

% running the module to create figure 4
% figure 4 is the baseline simulation
% figures are saved in folder: bubblewallfigs and contourfigs
run m_run_cavitation_baseline.m

% running the module to create figure 4
% figure 4 is the baseline simulation
% figures are saved in folder: bubblewallfigs and contourfigs
run m_run_cavitation_ML.m

% running the module to create figure 4
% figure 4 is the baseline simulation
% figures are saved in folder: bubblewallfigs and contourfigs
run m_run_cavitation_MN.m

% running the module to create figure 4
% figure 4 is the baseline simulation
% figures are saved in folder: bubblewallfigs and contourfigs
run m_run_cavitation_HN.m

% running the module to create figure 4
% figure 4 is the baseline simulation
% figures are saved in folder: bubblewallfigs and contourfigs

close all; clc; clear;
run m_run_cavitation_frequency.m
run m_run_cavitation_radius.m
run m_run_cavitation_pressure.m

close all; clc;