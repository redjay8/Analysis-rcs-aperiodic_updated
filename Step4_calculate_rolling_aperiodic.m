%% Step4_calculate_rolling_aperiodic.m
% Method 1 (60s windows with 5s steps):
% High temporal resolution, but multiple neural measures map to one PKG reading

% Method 2: Hybrid approach (120s steps, aligned to end at PKG times):
% Independent measurements with no window overlap
% Matches the sampling rate of both data streams

% Method 3: Averaging approach (5s steps, averaged per PKG interval):
% Leverages all data but reduces dimensionality

clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));
addpath(genpath('../code/fieldtrip')); % Ensure FieldTrip is in the path
ft_defaults;

%% 2. Load patient data and condition states
patient = 'Patient1';
side = 'L';
brain_regions = {'+2-0', '+3-1', '+9-8', '+11-10'}; 
brain_region_name = {'STN1', 'STN1', 'postcentral gyrus', 'precentral gyrus'}; 
chan = 2; % Target electrode channel for analysis
pkg_patient = 'Pat1';  % Corresponding PKG patient name (if different)

pathName = fullfile(base_path, patient, 'Device', 'Session');
combined_file = fullfile(pathName, 'combinedDataTable.mat');
state_file = fullfile(pathName, 'condition_states.mat');

fprintf('Loading neural data...\n');
load(combined_file);
load(state_file);

%% 3. Load PKG data to get actual PKG timestamps
fprintf('Loading PKG data to get actual timestamps...\n');

pkg_dir = fullfile(base_path, 'PKG_Data');
if ~exist(pkg_dir, 'dir')
    error('PKG data directory not found: %s', pkg_dir);
end

% Look for PKG files matching the patient name
pkg_files = dir(fullfile(pkg_dir, sprintf('%s_*.csv', pkg_patient)));
if isempty(pkg_files)
    error('No PKG CSV files found for patient %s. Please check if PKG data is available.', pkg_patient);
end

% Load and concatenate all PKG files for this patient
pkg_data = table();
for i = 1:length(pkg_files)
    fprintf('Loading PKG file %d/%d: %s\n', i, length(pkg_files), pkg_files(i).name);
    
    % Read CSV file
    file_path = fullfile(pkg_dir, pkg_files(i).name);
    pkg_file_data = readtable(file_path, 'TextType', 'string');
    
    % Append to combined data
    pkg_data = [pkg_data; pkg_file_data];
end

% Convert date-time to POSIX time for easier alignment
pkg_data.datetime = datetime(pkg_data.Date_Time);
pkg_data.time = posixtime(pkg_data.datetime);

% Sort by time for proper alignment
pkg_data = sortrows(pkg_data, 'time');

% Filter PKG times to those within the neural recording period plus a buffer
neural_start_time = combinedDataTable.DerivedTime(1);
neural_end_time = combinedDataTable.DerivedTime(end);
pkg_interval_sec = 120; % Standard 2-minute PKG interval
pkg_in_range = pkg_data.time >= neural_start_time & pkg_data.time <= (neural_end_time + pkg_interval_sec);
pkg_times = pkg_data.time(pkg_in_range);

fprintf('Found %d PKG timestamps within or just after neural recording range\n', length(pkg_times));

% If no PKG timestamps in range, use the estimated approach as fallback
if isempty(pkg_times)
    warning('No PKG timestamps found within neural recording period. Using estimated intervals.');
    pkg_start_time = neural_start_time + pkg_interval_sec; % First PKG timestamp after start
    pkg_end_time = neural_end_time + pkg_interval_sec;
    pkg_times = (pkg_start_time:pkg_interval_sec:pkg_end_time)';
    fprintf('Using %d estimated PKG timestamps\n', length(pkg_times));
end

%% 4. Set up parameters for analysis
fprintf('Setting up analysis parameters...\n');

% Channel selection - which LFP channel to analyze
channel_name = ['TD_key', num2str(chan-1)]; % e.g., TD_key0 for chan=1
if ~ismember(channel_name, combinedDataTable.Properties.VariableNames)
    error('Channel %s not found in combinedDataTable', channel_name);
end

% Time series of LFP data
time_vector = combinedDataTable.DerivedTime;
signal = combinedDataTable.(channel_name);

% Get sampling rate from timeDomainSettings
fs = timeDomainSettings.samplingRate(1);

% Common window and frequency parameters
window_size_sec = 60;     % 60-second window size, matching Wiest et al.
window_size = window_size_sec * fs;  % Convert to samples

% Define frequency ranges based on Wiest et al.
med_freq_range = [40 90];   % Range for medication analysis (avoid beta peaks)
stim_freq_range = [10 50];  % Range for DBS analysis (avoid spectral plateau)

%% METHODS 1,2,3 for aperiodic metrics

%% 5a. METHOD 1: Original method - 60s windows with 5s steps
fprintf('\n--- METHOD 1: Original method - 60s windows with 5s steps ---\n');

% Step size for original approach
step_size_sec = 5;        % 5-second step size
step_size = step_size_sec * fs;      % Convert to samples

% Pre-allocate arrays for results
num_windows = floor((length(signal) - window_size) / step_size) + 1;
aperiodic_times_m1 = zeros(num_windows, 1);
aperiodic_exponents_m1 = zeros(num_windows, 1);
aperiodic_offsets_m1 = zeros(num_windows, 1);
aperiodic_r_squared_m1 = zeros(num_windows, 1);
window_med_state_m1 = zeros(num_windows, 1);
window_stim_state_m1 = zeros(num_windows, 1);

% Process each window
for i = 1:num_windows
    % Extract window indices
    start_idx = (i-1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    
    % Skip if we reach the end of the data
    if end_idx > length(signal)
        break;
    end
    
    % Extract window data
    window_data = signal(start_idx:end_idx);
    window_time = time_vector(start_idx:end_idx);
    
    % Store the center time of this window
    aperiodic_times_m1(i) = window_time(round(length(window_time)/2));
    
    % Calculate window medication and stimulation states (most common state in window)
    window_med_state_m1(i) = mode(med_state(start_idx:end_idx));
    window_stim_state_m1(i) = mode(stim_state(start_idx:end_idx));
    
    % Skip windows with too many NaNs
    nan_idx = isnan(window_data);
    if sum(nan_idx) > window_size * 0.2 % Skip if >20% NaNs
        aperiodic_exponents_m1(i) = NaN;
        aperiodic_offsets_m1(i) = NaN;
        aperiodic_r_squared_m1(i) = NaN;
        continue;
    end
    
    window_data(nan_idx) = 0; % Replace NaNs with zeros for processing
    
    % Determine appropriate frequency range based on stimulation state
    if window_stim_state_m1(i) > 0.5
        freq_range = stim_freq_range;
    else
        freq_range = med_freq_range;
    end
    
    % Calculate aperiodic component and store results
    [exponent, offset, r_squared] = calculate_aperiodic(window_data, fs, freq_range, channel_name);
    
    aperiodic_exponents_m1(i) = exponent;
    aperiodic_offsets_m1(i) = offset;
    aperiodic_r_squared_m1(i) = r_squared;
    
    % Display progress
    if mod(i, 100) == 0
        fprintf('Method 1: Processed %d/%d windows (%.1f%%)\n', i, num_windows, 100*i/num_windows);
    end
end

% Remove windows with NaN exponents
valid_idx = ~isnan(aperiodic_exponents_m1);
aperiodic_times_m1 = aperiodic_times_m1(valid_idx);
aperiodic_exponents_m1 = aperiodic_exponents_m1(valid_idx);
aperiodic_offsets_m1 = aperiodic_offsets_m1(valid_idx);
aperiodic_r_squared_m1 = aperiodic_r_squared_m1(valid_idx);
window_med_state_m1 = window_med_state_m1(valid_idx);
window_stim_state_m1 = window_stim_state_m1(valid_idx);

%% 5b. METHOD 2: 60s windows ending exactly at PKG measurement times
fprintf('\n--- METHOD 2: 60s windows ending at PKG measurement times ---\n');

% Pre-allocate arrays for Method 2 results
num_pkg_intervals = length(pkg_times);
aperiodic_times_m2 = pkg_times; % End times match PKG times exactly
aperiodic_exponents_m2 = zeros(num_pkg_intervals, 1);
aperiodic_offsets_m2 = zeros(num_pkg_intervals, 1);
aperiodic_r_squared_m2 = zeros(num_pkg_intervals, 1);
window_med_state_m2 = zeros(num_pkg_intervals, 1);
window_stim_state_m2 = zeros(num_pkg_intervals, 1);

% Calculate and display the actual step size (typically 120s for PKG data)
if length(pkg_times) > 1
    pkg_intervals = diff(pkg_times);
    avg_pkg_interval = mean(pkg_intervals);
    fprintf('Method 2: Using PKG timestamps with average interval of %.1f seconds\n', avg_pkg_interval);
else
    fprintf('Method 2: Only one PKG timestamp available\n');
end

% Process each window ending at PKG timestamp
for i = 1:num_pkg_intervals
    % Calculate window boundaries: 60s window ending exactly at PKG time
    end_time = pkg_times(i);
    start_time = end_time - window_size_sec;
    
    % Find closest indices in the neural time vector
    [~, end_idx] = min(abs(time_vector - end_time));
    [~, start_idx] = min(abs(time_vector - start_time));
    
    % Skip if window is outside the neural data range
    if start_idx < 1 || end_idx > length(signal) || start_idx >= end_idx
        aperiodic_exponents_m2(i) = NaN;
        aperiodic_offsets_m2(i) = NaN;
        aperiodic_r_squared_m2(i) = NaN;
        continue;
    end
    
    % Extract window data
    window_data = signal(start_idx:end_idx);
    
    % Calculate window medication and stimulation states (most common state in window)
    window_med_state_m2(i) = mode(med_state(start_idx:end_idx));
    window_stim_state_m2(i) = mode(stim_state(start_idx:end_idx));
    
    % Skip windows with too many NaNs
    nan_idx = isnan(window_data);
    if sum(nan_idx) > length(window_data) * 0.2 % Skip if >20% NaNs
        aperiodic_exponents_m2(i) = NaN;
        aperiodic_offsets_m2(i) = NaN;
        aperiodic_r_squared_m2(i) = NaN;
        continue;
    end
    
    window_data(nan_idx) = 0; % Replace NaNs with zeros for processing
    
    % Determine appropriate frequency range based on stimulation state
    if window_stim_state_m2(i) > 0.5
        freq_range = stim_freq_range;
    else
        freq_range = med_freq_range;
    end
    
    % Calculate aperiodic component and store results
    [exponent, offset, r_squared] = calculate_aperiodic(window_data, fs, freq_range, channel_name);
    
    aperiodic_exponents_m2(i) = exponent;
    aperiodic_offsets_m2(i) = offset;
    aperiodic_r_squared_m2(i) = r_squared;
    
    % Display progress
    if mod(i, 10) == 0
        fprintf('Method 2: Processed %d/%d windows (%.1f%%)\n', i, num_pkg_intervals, 100*i/num_pkg_intervals);
    end
end

% Remove windows with NaN exponents
valid_idx = ~isnan(aperiodic_exponents_m2);
aperiodic_times_m2 = aperiodic_times_m2(valid_idx);
aperiodic_exponents_m2 = aperiodic_exponents_m2(valid_idx);
aperiodic_offsets_m2 = aperiodic_offsets_m2(valid_idx);
aperiodic_r_squared_m2 = aperiodic_r_squared_m2(valid_idx);
window_med_state_m2 = window_med_state_m2(valid_idx);
window_stim_state_m2 = window_stim_state_m2(valid_idx);

%% 5c. METHOD 3: Average exponents in 120s window immediately preceding PKG measurements
fprintf('\n--- METHOD 3: Average exponents in 120s window preceding PKG measurements ---\n');

% Pre-allocate arrays for Method 3 results
aperiodic_times_m3 = pkg_times; % End times match PKG times
aperiodic_exponents_m3 = zeros(num_pkg_intervals, 1);
aperiodic_offsets_m3 = zeros(num_pkg_intervals, 1);
aperiodic_r_squared_m3 = zeros(num_pkg_intervals, 1);
window_med_state_m3 = zeros(num_pkg_intervals, 1);
window_stim_state_m3 = zeros(num_pkg_intervals, 1);

% Calculate and display the actual step size (typically 120s for PKG data)
if length(pkg_times) > 1
    pkg_intervals = diff(pkg_times);
    avg_pkg_interval = mean(pkg_intervals);
    fprintf('Method 3: Using PKG timestamps with average interval of %.1f seconds\n', avg_pkg_interval);
else
    fprintf('Method 3: Only one PKG timestamp available\n');
end

% For each PKG timestamp, find all Method 1 windows that fall within the 120s window before it
for i = 1:num_pkg_intervals
    % Define the 120s window preceding the PKG timestamp
    pkg_time = pkg_times(i);
    window_start = pkg_time - 120; % 120s window (as PKG data is every 2 min)
    
    % Find Method 1 results within this window
    window_idx = aperiodic_times_m1 >= window_start & aperiodic_times_m1 <= pkg_time;
    
    % If no exponents in this window, use NaN
    if sum(window_idx) == 0
        aperiodic_exponents_m3(i) = NaN;
        aperiodic_offsets_m3(i) = NaN;
        aperiodic_r_squared_m3(i) = NaN;
        window_med_state_m3(i) = NaN;
        window_stim_state_m3(i) = NaN;
        continue;
    end
    
    % Average the values within the 120s window
    aperiodic_exponents_m3(i) = mean(aperiodic_exponents_m1(window_idx), 'omitnan');
    aperiodic_offsets_m3(i) = mean(aperiodic_offsets_m1(window_idx), 'omitnan');
    aperiodic_r_squared_m3(i) = mean(aperiodic_r_squared_m1(window_idx), 'omitnan');
    
    % For states, take the mode (most common state) within the window
    window_med_state_m3(i) = mode(window_med_state_m1(window_idx));
    window_stim_state_m3(i) = mode(window_stim_state_m1(window_idx));
    
    % Display progress
    if mod(i, 10) == 0
        fprintf('Method 3: Processed %d/%d intervals (%.1f%%)\n', i, num_pkg_intervals, 100*i/num_pkg_intervals);
    end
end

% Remove windows with NaN exponents
valid_idx = ~isnan(aperiodic_exponents_m3);
aperiodic_times_m3 = aperiodic_times_m3(valid_idx);
aperiodic_exponents_m3 = aperiodic_exponents_m3(valid_idx);
aperiodic_offsets_m3 = aperiodic_offsets_m3(valid_idx);
aperiodic_r_squared_m3 = aperiodic_r_squared_m3(valid_idx);
window_med_state_m3 = window_med_state_m3(valid_idx);
window_stim_state_m3 = window_stim_state_m3(valid_idx);

%% 6. Save results from all methods
fprintf('\nSaving aperiodic exponent time series from all methods...\n');

% Organize results in tables
aperiodic_results_m1 = table();
aperiodic_results_m1.time = aperiodic_times_m1;
aperiodic_results_m1.exponent = aperiodic_exponents_m1;
aperiodic_results_m1.offset = aperiodic_offsets_m1;
aperiodic_results_m1.r_squared = aperiodic_r_squared_m1;
aperiodic_results_m1.med_state = window_med_state_m1;
aperiodic_results_m1.stim_state = window_stim_state_m1;

aperiodic_results_m2 = table();
aperiodic_results_m2.time = aperiodic_times_m2;
aperiodic_results_m2.exponent = aperiodic_exponents_m2;
aperiodic_results_m2.offset = aperiodic_offsets_m2;
aperiodic_results_m2.r_squared = aperiodic_r_squared_m2;
aperiodic_results_m2.med_state = window_med_state_m2;
aperiodic_results_m2.stim_state = window_stim_state_m2;

aperiodic_results_m3 = table();
aperiodic_results_m3.time = aperiodic_times_m3;
aperiodic_results_m3.exponent = aperiodic_exponents_m3;
aperiodic_results_m3.offset = aperiodic_offsets_m3;
aperiodic_results_m3.r_squared = aperiodic_r_squared_m3;
aperiodic_results_m3.med_state = window_med_state_m3;
aperiodic_results_m3.stim_state = window_stim_state_m3;

% Add metadata
metadata = struct();
metadata.patient = patient;
metadata.side = side;
metadata.brain_region = brain_regions{chan};
metadata.brain_region_name = brain_region_name{chan};
metadata.window_size_sec = window_size_sec;
metadata.med_freq_range = med_freq_range;
metadata.stim_freq_range = stim_freq_range;
metadata.fooof_implementation = 'FieldTrip';
metadata.pkg_interval_sec = pkg_interval_sec;
metadata.pkg_patient = pkg_patient;

% Add method-specific metadata
metadata.m1_description = 'Method 1: 60s windows with 5s steps';
metadata.m2_description = 'Method 2: 60s windows ending exactly at PKG times';
metadata.m3_description = 'Method 3: Average exponents in 120s window preceding PKG measurements';
metadata.m1_step_size_sec = step_size_sec;

% Save to combined file
results_file = fullfile(pathName, sprintf('aperiodic_timeseries_all_%s_%s.mat', brain_regions{chan}, side));
save(results_file, 'aperiodic_results_m1', 'aperiodic_results_m2', 'aperiodic_results_m3', 'metadata', 'pkg_times');
fprintf('All methods aperiodic time series saved to: %s\n', results_file);

% Save individual method files for backward compatibility
results_file_m1 = fullfile(pathName, sprintf('aperiodic_timeseries_%s_%s.mat', brain_regions{chan}, side));
aperiodic_results = aperiodic_results_m1; % For backward compatibility
save(results_file_m1, 'aperiodic_results', 'metadata');

results_file_m2 = fullfile(pathName, sprintf('aperiodic_timeseries_m2_%s_%s.mat', brain_regions{chan}, side));
aperiodic_results = aperiodic_results_m2;
save(results_file_m2, 'aperiodic_results', 'metadata');

results_file_m3 = fullfile(pathName, sprintf('aperiodic_timeseries_m3_%s_%s.mat', brain_regions{chan}, side));
aperiodic_results = aperiodic_results_m3;
save(results_file_m3, 'aperiodic_results', 'metadata');

%% 7. Visualize and compare the three methods
fprintf('Visualizing and comparing all methods...\n');

% Create a figure to compare all methods
figure('Position', [100 100 1200 800]);

% Plot aperiodic exponents time series comparison
ax1 = subplot(3,1,1);
plot(datetime(aperiodic_times_m1, 'ConvertFrom', 'posixtime'), aperiodic_exponents_m1, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
hold on;
plot(datetime(aperiodic_times_m2, 'ConvertFrom', 'posixtime'), aperiodic_exponents_m2, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
plot(datetime(aperiodic_times_m3, 'ConvertFrom', 'posixtime'), aperiodic_exponents_m3, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg before PKG)');
title('Aperiodic Exponent Time Series Comparison');
ylabel('Exponent');
legend('Location', 'best');
grid on;

% Plot aperiodic offset time series comparison
ax2 = subplot(3,1,2);
plot(datetime(aperiodic_times_m1, 'ConvertFrom', 'posixtime'), aperiodic_offsets_m1, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
hold on;
plot(datetime(aperiodic_times_m2, 'ConvertFrom', 'posixtime'), aperiodic_offsets_m2, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
plot(datetime(aperiodic_times_m3, 'ConvertFrom', 'posixtime'), aperiodic_offsets_m3, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg before PKG)');
title('Aperiodic Offset Time Series Comparison');
ylabel('Offset');
grid on;

% Plot goodness of fit comparison
ax3 = subplot(3,1,3);
plot(datetime(aperiodic_times_m1, 'ConvertFrom', 'posixtime'), aperiodic_r_squared_m1, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
hold on;
plot(datetime(aperiodic_times_m2, 'ConvertFrom', 'posixtime'), aperiodic_r_squared_m2, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
plot(datetime(aperiodic_times_m3, 'ConvertFrom', 'posixtime'), aperiodic_r_squared_m3, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg before PKG)');
title('Goodness of Fit (R^2) Comparison');
ylabel('R^2');
xlabel('Time');
grid on;
ylim([0 1]);

% Link axes for synchronized zooming
linkaxes([ax1, ax2, ax3], 'x');

% Save figure
saveas(gcf, fullfile(pathName, sprintf('aperiodic_timeseries_comparison_%s_%s.fig', brain_regions{chan}, side)));
saveas(gcf, fullfile(pathName, sprintf('aperiodic_timeseries_comparison_%s_%s.png', brain_regions{chan}, side)));

% Calculate and display statistics for each method
fprintf('\nMethod statistics:\n');
fprintf('Method 1: %d measurements, mean exponent = %.4f, SD = %.4f\n', ...
    length(aperiodic_exponents_m1), mean(aperiodic_exponents_m1), std(aperiodic_exponents_m1));
fprintf('Method 2: %d measurements, mean exponent = %.4f, SD = %.4f\n', ...
    length(aperiodic_exponents_m2), mean(aperiodic_exponents_m2), std(aperiodic_exponents_m2));
fprintf('Method 3: %d measurements, mean exponent = %.4f, SD = %.4f\n', ...
    length(aperiodic_exponents_m3), mean(aperiodic_exponents_m3), std(aperiodic_exponents_m3));

fprintf('Step 4 processing complete with all three methods.\n');

%% Helper function to calculate aperiodic component using FieldTrip

function [exponent, offset, r_squared] = calculate_aperiodic(window_data, fs, freq_range, channel_name)
    % Prepare data for FieldTrip processing
    ft_data = [];
    ft_data.trial{1} = window_data';
    ft_data.time{1} = (0:length(window_data)-1) / fs;
    ft_data.label = {channel_name};
    ft_data.fsample = fs;
    
    % Calculate aperiodic component using FieldTrip's FOOOF implementation
    cfg = [];
    cfg.foilim = freq_range;  % Use appropriate frequency range
    cfg.pad = 'maxperlen';
    cfg.tapsmofrq = 2;
    cfg.method = 'mtmfft';
    cfg.output = 'fooof_aperiodic';
    % FOOOF specific settings
    cfg.fooof.peak_width_limits = [2 8];
    cfg.fooof.max_n_peaks = 6;
    cfg.fooof.min_peak_height = 0.1;
    cfg.fooof.peak_threshold = 2.0;
    cfg.fooof.aperiodic_mode = 'fixed';  % Use fixed mode for consistent linear fitting
    
    try
        aperiodic = ft_freqanalysis(cfg, ft_data);
        
        % Extract aperiodic parameters from FieldTrip's output
        if isfield(aperiodic, 'aperiodic')
            % Extract aperiodic exponent (slope) and offset
            exponent = aperiodic.aperiodic.exponent;
            offset = aperiodic.aperiodic.offset;
            
            % Calculate goodness of fit (R²)
            if isfield(aperiodic.aperiodic, 'r_squared')
                r_squared = aperiodic.aperiodic.r_squared;
            else
                % If not available directly, set to NaN
                r_squared = NaN;
            end
        else
            % If FOOOF processing failed but didn't raise an error
            exponent = NaN;
            offset = NaN;
            r_squared = NaN;
        end
    catch
        % If an error occurred during FOOOF processing
        exponent = NaN;
        offset = NaN;
        r_squared = NaN;
    end
end