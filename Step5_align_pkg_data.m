%% Step5_align_pkg_data.m
% Import PKG data, sort labeled states, and align with aperiodic exponents
% Based on Code Ocean
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));

%% 2. Load aperiodic data
patient = 'Patient1';  % RC+S recording patient name
pkg_patient = 'Pat1';  % Corresponding PKG patient name (if different)
side = 'L';
brain_regions = {'+2-0', '+3-1', '+9-8', '+11-10'}; 
brain_region_name = {'STN1', 'STN1', 'postcentral gyrus', 'precentral gyrus'}; 
chan = 2;

% Load aperiodic exponent time series from Step 4
pathName = fullfile(base_path, patient, 'Device', 'Session');
aperiodic_file = fullfile(pathName, sprintf('aperiodic_timeseries_%s_%s.mat', brain_regions{chan}, side));
fprintf('Loading aperiodic results from: %s\n', aperiodic_file);
load(aperiodic_file);

%% 3. Load PKG data using your existing loading process
fprintf('Loading PKG data...\n');

% Map side abbreviation to full name for PKG filename
if strcmpi(side, 'L')
    hemi_full = 'Left';
    ophemi_full = 'Right';
else
    hemi_full = 'Right';
    ophemi_full = 'Left';
end

pkg_file = fullfile(base_path, sprintf('%s_%sBrain_%sBody.mat', pkg_patient, hemi_full, ophemi_full));
fprintf('Loading PKG data from: %s\n', pkg_file);

try
    % Load the raw PKG data file
    load(pkg_file);
    
    % Set up PKG structure
    cS = pkg_patient;
    cH = side;
    cCo = [cS cH];
    
    % Check if the structure needs initialization
    if ~isfield(pkg, cCo) || ~istable(pkg.(cCo))
        pkg.(cCo) = table;
        
        % Set variable names
        pkg.(cCo).Properties.VariableNames = {'ldState', 'DKS', 'BKS', 'TS', 'cond', 'daynum'};
    end
    
    % Process state values
    pkg.(cCo).ldState(pkg.(cCo).ldState == 3) = 0;
    pkg.(cCo).ldState(pkg.(cCo).ldState == 4) = 1;
    pkg.(cCo).ldState(pkg.(cCo).ldState == 5) = 2;
    pkg.(cCo).ldState(pkg.(cCo).ldState == 8) = 1;
    
    pkg_data_available = true;
    fprintf('Successfully loaded PKG data\n');
    
end

%% 4. Process PKG data
% To-do, pending data structure
