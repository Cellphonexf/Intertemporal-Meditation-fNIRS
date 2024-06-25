%% fNIRS data visualization
%% Programmed by Feng Xiao (2024.6.25)
clear all,
clc,
%% Load fNIRS_visualization.xlsx
[~, sheetNames] = xlsfinfo('fNIRS_visualization.xlsx');

for i = 1:numel(sheetNames)
    [data, ~, ~] = xlsread('fNIRS_visualization.xlsx', sheetNames{i});
    assignin('base', sheetNames{i}, data);
end
df = 51;
%% Meditation contrast (session 1)
HbO2 = contrast_ses1;
Hb = ref; %because EasyTopo is required to input both the HbO2 and Hb for visualization, so we used the sham-data (i.e., all data points were -1) as reference
save('input_contrast1_data.mat', 'mni', 'HbO2', 'Hb', 'df');
%% Meditation contrast (session 1 & 2)
HbO2 = contrast_ses12;
Hb = ref; %because EasyTopo is required to input both the HbO2 and Hb for visualization, so we used the sham-data (i.e., all data points were -1) as reference
save('input_contrast12_data.mat', 'mni', 'HbO2', 'Hb', 'df');
%% Meditation contrast (session 1 & 2 & 3)
HbO2 = contrast_ses123;
Hb = ref; %because EasyTopo is required to input both the HbO2 and Hb for visualization, so we used the sham-data (i.e., all data points were -1) as reference
save('input_contrast123_data.mat', 'mni', 'HbO2', 'Hb', 'df');