clear all; close all; clc;

%reading the file name
filename = strcat('PolyAcry_S2_15_48_42.mat');
% loading the file
load(filename);
% converting the structure to cell, you can not do this from the internal
% structure object
Cell = struct2cell(BubblePlotSet);
% First object in the structure is the radius
% Accessing the radius data
R = Cell{1}';
time = Cell{14}';
% Removing NaNs to now write the function
idx = ~isnan(R);
time = time(idx);
R = R(idx);
% From here you can print out the CSV file

