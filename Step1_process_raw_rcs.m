%% Step1_process_raw_rcs.m
% Script to process raw RC+S JSON files using ProcessRCS
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));

%% 2. Process raw RC+S data with ProcessRCS
% Full path to RC+S Device folder with raw JSON files
pathName = fullfile(base_path, 'Patient1', 'Device', 'Session');

% Process flag options:
% 1: Process and save (overwrite if processed file already exists)
% 2: Process and do not save
% 3: If processed file already exists, then load. If not, process and save
% 4: If processed file already exists, then load. If not, process but do not save
processFlag = 3; 
shortGaps_systemTick = 0;

% Process the data
fprintf('Processing RC+S data from: %s\n', pathName);
[unifiedDerivedTimes, timeDomainData, timeDomainData_onlyTimeVariables, timeDomain_timeVariableNames,...
 AccelData, AccelData_onlyTimeVariables, Accel_timeVariableNames,...
 PowerData, PowerData_onlyTimeVariables, Power_timeVariableNames,...
 FFTData, FFTData_onlyTimeVariables, FFT_timeVariableNames,...
 AdaptiveData, AdaptiveData_onlyTimeVariables, Adaptive_timeVariableNames,...
 timeDomainSettings, powerSettings, fftSettings, eventLogTable, metaData,...
 stimSettingsOut, stimMetaData, stimLogSettings,...
 DetectorSettings, AdaptiveStimSettings, AdaptiveEmbeddedRuns_StimSettings] = ...
    ProcessRCS(pathName, processFlag, shortGaps_systemTick);

fprintf('Step 1 processing complete.\n');