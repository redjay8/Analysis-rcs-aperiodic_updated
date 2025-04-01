%% Step2_create_combined_table.m
% Create combined data table from processed RC+S data
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));

%% 2. Load processed data (the output of ProcessRCS)
pathName = fullfile(base_path, 'Patient1', 'Device', 'Session');
load_file = fullfile(pathName, 'AllDataTables.mat');
fprintf('Loading data from: %s\n', load_file);
load(load_file);

%% 3. Create combined data table following the repository's recommended approach
fprintf('Creating combined data table...\n');
dataStreams = {timeDomainData, AccelData, PowerData, FFTData, AdaptiveData};
[combinedDataTable] = createCombinedTable(dataStreams, unifiedDerivedTimes, metaData);

%% 4. Save combined data table
save_file = fullfile(pathName, 'combinedDataTable.mat');
save(save_file, 'combinedDataTable', 'timeDomainSettings', 'eventLogTable', 'stimSettingsOut', 'AdaptiveStimSettings');
fprintf('Combined data table saved to: %s\n', save_file);

fprintf('Step 2 processing complete.\n');