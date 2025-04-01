%% Step4_calculate_rolling_aperiodic.m
% Calculate aperiodic exponents using a rolling window approach
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));
addpath(genpath('../code/fieldtrip')); % Ensure FieldTrip is in the path
ft_defaults;

%% 2. Load patient data and condition states [TENTATIVE - Code Ocean defaults]
patient = 'Patient1';
side = 'L';
brain_regions = {'+2-0', '+3-1', '+9-8', '+11-10'}; % Example channel pairs 
brain_region_name = {'STN1', 'STN1', 'postcentral gyrus', 'precentral gyrus'}; 
chan = 2; % Target electrode channel for analysis

pathName = fullfile(base_path, patient, 'Device', 'Session');
combined_file = fullfile(pathName, 'combinedDataTable.mat');
state_file = fullfile(pathName, 'condition_states.mat');

fprintf('Loading data...\n');
load(combined_file);
load(state_file);

%% 3. Set up parameters for rolling window analysis
fprintf('Setting up rolling window parameters...\n');

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

% Window parameters
window_size_sec = 60;     % 60-second window size, matching Wiest et al.
window_size = window_size_sec * fs;  % Convert to samples
step_size_sec = 5;        % 5-second step size
step_size = step_size_sec * fs;      % Convert to samples

% Define frequency ranges based on Wiest et al.
med_freq_range = [40 90];   % Range for medication analysis (avoid beta peaks)
stim_freq_range = [10 50];  % Range for DBS analysis (avoid spectral plateau)

%% 4. Calculate aperiodic exponents using FieldTrip's FOOOF implementation
fprintf('Calculating aperiodic exponents using FieldTrip''s FOOOF implementation...\n');

% Pre-allocate arrays for results
num_windows = floor((length(signal) - window_size) / step_size) + 1;
aperiodic_times = zeros(num_windows, 1);
aperiodic_exponents = zeros(num_windows, 1);
aperiodic_offsets = zeros(num_windows, 1);
aperiodic_r_squared = zeros(num_windows, 1);
window_med_state = zeros(num_windows, 1);
window_stim_state = zeros(num_windows, 1);

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
    aperiodic_times(i) = window_time(round(length(window_time)/2));
    
    % Calculate window medication and stimulation states (most common state in window)
    window_med_state(i) = mode(med_state(start_idx:end_idx));
    window_stim_state(i) = mode(stim_state(start_idx:end_idx));
    
    % Skip windows with too many NaNs
    nan_idx = isnan(window_data);
    if sum(nan_idx) > window_size * 0.2 % Skip if >20% NaNs
        aperiodic_exponents(i) = NaN;
        aperiodic_offsets(i) = NaN;
        aperiodic_r_squared(i) = NaN;
        continue;
    end
    
    window_data(nan_idx) = 0; % Replace NaNs with zeros for processing
    
    % Determine appropriate frequency range based on stimulation state
    if window_stim_state(i) > 0.5
        freq_range = stim_freq_range;
    else
        freq_range = med_freq_range;
    end
    
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
            aperiodic_exponents(i) = aperiodic.aperiodic.exponent;
            aperiodic_offsets(i) = aperiodic.aperiodic.offset;
            
            % Calculate goodness of fit (R²)
            if isfield(aperiodic.aperiodic, 'r_squared')
                aperiodic_r_squared(i) = aperiodic.aperiodic.r_squared;
            else
                % If not available directly, we could calculate it ourselves
                % But for now, just set to NaN
                aperiodic_r_squared(i) = NaN;
            end
        else
            % If FOOOF processing failed but didn't raise an error
            aperiodic_exponents(i) = NaN;
            aperiodic_offsets(i) = NaN;
            aperiodic_r_squared(i) = NaN;
        end
    catch ME
        % If an error occurred during FOOOF processing
        warning('FOOOF processing failed for window %d: %s', i, ME.message);
        aperiodic_exponents(i) = NaN;
        aperiodic_offsets(i) = NaN;
        aperiodic_r_squared(i) = NaN;
    end
    
    % Display progress
    if mod(i, 10) == 0
        fprintf('Processed %d/%d windows (%.1f%%)\n', i, num_windows, 100*i/num_windows);
    end
end

% Remove windows with NaN exponents
valid_idx = ~isnan(aperiodic_exponents);
aperiodic_times = aperiodic_times(valid_idx);
aperiodic_exponents = aperiodic_exponents(valid_idx);
aperiodic_offsets = aperiodic_offsets(valid_idx);
aperiodic_r_squared = aperiodic_r_squared(valid_idx);
window_med_state = window_med_state(valid_idx);
window_stim_state = window_stim_state(valid_idx);

%% 5. Save results
fprintf('Saving aperiodic exponent time series...\n');

% Organize results in a table
aperiodic_results = table();
aperiodic_results.time = aperiodic_times;
aperiodic_results.exponent = aperiodic_exponents;
aperiodic_results.offset = aperiodic_offsets;
aperiodic_results.r_squared = aperiodic_r_squared;
aperiodic_results.med_state = window_med_state;
aperiodic_results.stim_state = window_stim_state;

% Add metadata
metadata = struct();
metadata.patient = patient;
metadata.side = side;
metadata.brain_region = brain_regions{chan};
metadata.brain_region_name = brain_region_name{chan};
metadata.window_size_sec = window_size_sec;
metadata.step_size_sec = step_size_sec;
metadata.med_freq_range = med_freq_range;
metadata.stim_freq_range = stim_freq_range;
metadata.fooof_implementation = 'FieldTrip';

% Save to file
results_file = fullfile(pathName, sprintf('aperiodic_timeseries_%s_%s.mat', brain_regions{chan}, side));
save(results_file, 'aperiodic_results', 'metadata');
fprintf('Aperiodic time series saved to: %s\n', results_file);

%% 6. Visualize results
figure('Position', [100 100 1200 800]);

% Plot aperiodic exponents time series
subplot(3,1,1);
plot(datetime(aperiodic_times, 'ConvertFrom', 'posixtime'), aperiodic_exponents, 'LineWidth', 1.5);
hold on;
% Highlight medication ON periods
med_on_idx = find(window_med_state == 1);
if ~isempty(med_on_idx)
    med_on_times = datetime(aperiodic_times(med_on_idx), 'ConvertFrom', 'posixtime');
    scatter(med_on_times, aperiodic_exponents(med_on_idx), 'g', 'filled', 'MarkerFaceAlpha', 0.3);
end
% Highlight stimulation ON periods
stim_on_idx = find(window_stim_state == 1);
if ~isempty(stim_on_idx)
    stim_on_times = datetime(aperiodic_times(stim_on_idx), 'ConvertFrom', 'posixtime');
    scatter(stim_on_times, aperiodic_exponents(stim_on_idx), 'r', 'filled', 'MarkerFaceAlpha', 0.3);
end
title('Aperiodic Exponent Time Series');
ylabel('Exponent');
grid on;

% Plot aperiodic offset time series
subplot(3,1,2);
plot(datetime(aperiodic_times, 'ConvertFrom', 'posixtime'), aperiodic_offsets, 'LineWidth', 1.5);
title('Aperiodic Offset Time Series');
ylabel('Offset');
grid on;

% Plot goodness of fit
subplot(3,1,3);
plot(datetime(aperiodic_times, 'ConvertFrom', 'posixtime'), aperiodic_r_squared, 'LineWidth', 1.5);
title('Goodness of Fit (R^2)');
ylabel('R^2');
xlabel('Time');
grid on;
ylim([0 1]);

% Save figure
saveas(gcf, fullfile(pathName, sprintf('aperiodic_timeseries_plot_%s_%s.fig', brain_regions{chan}, side)));
saveas(gcf, fullfile(pathName, sprintf('aperiodic_timeseries_plot_%s_%s.png', brain_regions{chan}, side)));

fprintf('Step 4 processing complete.\n');