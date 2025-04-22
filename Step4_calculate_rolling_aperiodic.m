%% Step4_calculate_rolling_aperiodic.m
%
% Calculates aperiodic features (exponent, offset) from LFP data using
% rolling windows and direct FOOOF Python implementation. Compares three
% different windowing strategies, potentially aligned with PKG timestamps.
%
% Method 1: Fixed 60s windows, 5s steps. High temporal resolution.
% Method 2: Fixed 60s windows, ending exactly at PKG timestamps. Aligned.
% Method 3: Average Method 1 results in 120s window preceding PKG times.
%
% **Assumptions:**
% 1. Steps 1-3 completed successfully.
% 2. Specific combinedDataTable & condition_states files exist in processed_data.
% 3. FieldTrip toolbox is available at the specified path.
% 4. FOOOF Python package is installed and accessible from MATLAB.
% 5. PKG data (if used) is available at the specified path and has a
%    'Date_Time' column readable by datetime().
% 6. The specified channel exists in the combinedDataTable.
% 7. **Sampling Rate Handling:** Attempts to load from timeDomainSettings,
%    but falls back to a hardcoded value (250 Hz for this dataset) if needed.
%
% **Output:**
% - Saves '.mat' files in 'processed_data' containing aperiodic results
%   tables and metadata for each method.
% - Saves comparison figures (.fig, .png) in 'processed_data'.

clear;
close all;
clc;

%% --- Configuration ---

% Define project and toolbox paths
project_base_path = '.\Code preparation';
toolbox_path = '.\Analysis-rcs-data-master';
fieldtrip_path = '.\fieldtrip-master'; % <--- *** UPDATE THIS PATH ***
pkg_data_folder = '.\PKG_Data'; % <--- *** UPDATE THIS PATH ***

% Define the folder for processed data output/input
processed_data_folder = fullfile(project_base_path, 'processed_data');

% --- Define which dataset to process (MUST MATCH PREVIOUS STEPS) ---
hemi_to_process = 'L';
session_id_to_process = 'Session1691533893108';
device_id_to_process = 'DeviceNPC700545H';

% --- Analysis Parameters ---
% Define channel mapping (ensure 'chan' is a valid index for these lists)
brain_regions = {'+2-0', '+3-1', '+9-8', '+11-10'}; % Example bipolar channels from DeviceSettings
brain_region_names = {'STN_L_02_00', 'STN_L_03_01', 'M1_L_09_08', 'M1_L_11_10'}; % Corresponding descriptive names
chan = 2; % Which channel index from the lists above to analyze (e.g., 2 for '+3-1')

% PKG related settings
pkg_patient_id = 'Pat1';  % Identifier used in PKG CSV filenames (e.g., 'Pat1_*.csv')
pkg_interval_sec = 120;   % Standard 2-minute PKG interval

% FOOOF / Windowing Parameters
window_size_sec = 60;     % Window size for spectral analysis
step_size_sec_m1 = 5;     % Step size for Method 1
med_freq_range = [40 90]; % Freq range for FOOOF when OFF/low stim
stim_freq_range = [10 50];% Freq range for FOOOF when ON/high stim (adjust based on stim artifact/beta changes)

% --- Sampling Rate Fallback ---
known_td_sampling_rate_hz = 250; % Hardcoded based on DeviceSettings.json analysis for this dataset

%% --- 1. Setup paths ---
fprintf('Starting Step 4: Calculate Rolling Aperiodic Metrics...\n');

% Add toolboxes to path
if ~contains(path, toolbox_path)
    addpath(genpath(toolbox_path));
    fprintf('Analysis-rcs-data toolbox added from: %s\n', toolbox_path);
end
if exist(fieldtrip_path, 'dir') && ~contains(path, fieldtrip_path)
    addpath(fieldtrip_path);
    ft_defaults; % Initialize FieldTrip
    fprintf('FieldTrip toolbox added from: %s\n', fieldtrip_path);
elseif ~exist(fieldtrip_path, 'dir')
     error('FieldTrip path not found: %s', fieldtrip_path);
else
    fprintf('FieldTrip toolbox already on path.\n');
    % Re-run ft_defaults just in case paths got shadowed
    ft_defaults; % Optional: uncomment if you suspect path issues
end

% Check if Python FOOOF module is accessible
try
    fooof_module = py.importlib.import_module('fooof');
    fprintf('FOOOF Python module found and accessible.\n');
catch ME
    warning('FOOOF Python module not found or not accessible: %s\nTrying to install it...', ME.message);
    try
        % Try to install FOOOF via pip
        cmd = sprintf('"%s" -m pip install fooof', pyversion);
        [status, cmdout] = system(cmd);
        if status == 0
            fprintf('FOOOF package installed successfully.\n');
            % Try to import again
            fooof_module = py.importlib.import_module('fooof');
            fprintf('FOOOF Python module now accessible.\n');
        else
            error('Failed to install FOOOF package. Error: %s\nPlease install manually using pip.', cmdout);
        end
    catch ME2
        error('Failed to access or install FOOOF Python module. Error: %s\nPlease install manually using pip.', ME2.message);
    end
end

%% --- 2. Load processed neural data and condition states ---

% Construct filenames based on configuration
combined_filename = sprintf('combinedDataTable_%s_%s_%s.mat', ...
                         hemi_to_process, session_id_to_process, device_id_to_process);
states_filename = sprintf('condition_states_%s_%s_%s.mat', ...
                          hemi_to_process, session_id_to_process, device_id_to_process);

combined_file_path = fullfile(processed_data_folder, combined_filename);
states_file_path = fullfile(processed_data_folder, states_filename);

% Load combined data table file
if ~isfile(combined_file_path)
    error('Combined data file not found: %s\nPlease ensure Step 2 ran successfully.', combined_file_path);
end
fprintf('Loading combined neural data from: %s\n', combined_file_path);
loaded_combined = load(combined_file_path); % Load into a struct

% Load condition states file
if ~isfile(states_file_path)
    error('Condition states file not found: %s\nPlease ensure Step 3 ran successfully.', states_file_path);
end
fprintf('Loading condition states from: %s\n', states_file_path);
loaded_states = load(states_file_path); % Load into a struct

% --- Verify required variables are loaded ---
% Note: We handle timeDomainSettings separately below due to potential loading issues
required_combined_vars = {'combinedDataTable'}; % Removed timeDomainSettings from strict check here
for v = 1:length(required_combined_vars)
    if ~isfield(loaded_combined, required_combined_vars{v})
        error('Variable "%s" not found in loaded file: %s', required_combined_vars{v}, combined_file_path);
    end
    % Assign to workspace
    eval([required_combined_vars{v} ' = loaded_combined.(required_combined_vars{v});']);
end

required_states_vars = {'med_state', 'stim_state'};
for v = 1:length(required_states_vars)
    if ~isfield(loaded_states, required_states_vars{v})
        error('Variable "%s" not found in loaded file: %s', required_states_vars{v}, states_file_path);
    end
     % Assign to workspace
    eval([required_states_vars{v} ' = loaded_states.(required_states_vars{v});']);
end
fprintf('Data loaded successfully.\n');

%% --- 3. Load PKG data (if available) ---
fprintf('Attempting to load PKG data...\n');
pkg_times = []; % Initialize as empty

if ~exist(pkg_data_folder, 'dir')
    warning('PKG data directory not found: %s. Cannot load PKG timestamps.', pkg_data_folder);
else
    % Look for PKG files matching the patient identifier
    pkg_file_pattern = sprintf('%s_*.csv', pkg_patient_id);
    pkg_files = dir(fullfile(pkg_data_folder, pkg_file_pattern));

    if isempty(pkg_files)
        warning('No PKG CSV files found matching pattern "%s" in %s.', pkg_file_pattern, pkg_data_folder);
    else
        % Load and concatenate all PKG files for this patient
        pkg_data_all = table();
        fprintf('Found %d PKG files. Loading...\n', length(pkg_files));
        for i = 1:length(pkg_files)
            file_path = fullfile(pkg_data_folder, pkg_files(i).name);
            fprintf('  Loading PKG file %d/%d: %s\n', i, length(pkg_files), pkg_files(i).name);
            try
                % Attempt to read CSV, be robust to potential read errors
                opts = detectImportOptions(file_path, 'TextType', 'string', 'ReadVariableNames', true);
                % Ensure 'Date_Time' column is read correctly
                var_idx = find(strcmp(opts.VariableNames, 'Date_Time'), 1);
                if isempty(var_idx)
                   warning('Skipping file %s: "Date_Time" column not found.', pkg_files(i).name);
                   continue;
                end
                opts = setvartype(opts, 'Date_Time', 'string'); % Ensure it's read as string first
                pkg_file_data = readtable(file_path, opts);

                % Check if 'Date_Time' column exists after reading
                 if ~ismember('Date_Time', pkg_file_data.Properties.VariableNames)
                     warning('Skipping file %s: "Date_Time" column missing after read.', pkg_files(i).name);
                     continue;
                 end

                % Convert date-time string to datetime objects, handle potential format issues
                try
                    pkg_file_data.datetime = datetime(pkg_file_data.Date_Time); % Try default format
                catch ME_datetime
                    warning('Could not parse Date_Time in file %s using default format. Error: %s. Trying common alternatives...', pkg_files(i).name, ME_datetime.message);
                    % Add alternative formats if needed, e.g.:
                    % pkg_file_data.datetime = datetime(pkg_file_data.Date_Time, 'InputFormat', 'MM/dd/yyyy HH:mm:ss');
                    % If conversion fails repeatedly, skip file or handle error
                    warning('Skipping file %s due to persistent Date_Time parsing issues.', pkg_files(i).name);
                    continue; % Skip this file if conversion fails
                end

                % Convert to POSIX time (seconds)
                pkg_file_data.time = posixtime(pkg_file_data.datetime);

                % Append to combined data (only necessary columns)
                 pkg_data_all = [pkg_data_all; pkg_file_data(:, {'time', 'datetime'})]; % Keep only needed cols

            catch ME_read
                warning('Error reading PKG file %s: %s. Skipping this file.', pkg_files(i).name, ME_read.message);
            end
        end % End loop through PKG files

        if ~isempty(pkg_data_all)
            % Sort by time
            pkg_data_all = sortrows(pkg_data_all, 'time');
            pkg_data_all = unique(pkg_data_all,'rows'); % Remove duplicate entries if files overlap

            % Filter PKG times to those relevant to the neural recording period
            % Convert neural times from ms (DerivedTime) to seconds for comparison
            neural_start_time_sec = combinedDataTable.DerivedTime(1) / 1000;
            neural_end_time_sec = combinedDataTable.DerivedTime(end) / 1000;
            % Include PKG times that end *during* or shortly *after* the neural recording
            pkg_in_range_idx = pkg_data_all.time >= neural_start_time_sec & pkg_data_all.time <= (neural_end_time_sec + pkg_interval_sec);
            pkg_times = pkg_data_all.time(pkg_in_range_idx);

            fprintf('Found %d unique PKG timestamps within or just after neural recording range (%.1f min duration).\n', ...
                    length(pkg_times), (neural_end_time_sec - neural_start_time_sec)/60);
        else
             warning('No valid PKG data loaded after checking all files.');
        end
    end % End check if PKG files exist
end % End check if PKG folder exists

% Fallback if no valid PKG times were loaded
if isempty(pkg_times)
    warning('No valid PKG timestamps loaded or found in range. Methods 2 and 3 will not be aligned to actual PKG data.');
    % Optional: Could generate estimated times here if needed for testing M2/M3 structure
    % pkg_start_time_sec = combinedDataTable.DerivedTime(1)/1000 + pkg_interval_sec;
    % pkg_end_time_sec = combinedDataTable.DerivedTime(end)/1000 + pkg_interval_sec;
    % pkg_times = (pkg_start_time_sec:pkg_interval_sec:pkg_end_time_sec)';
    % fprintf('Using %d estimated PKG timestamps for Methods 2 & 3 structure.\n', length(pkg_times));
    pkg_times = []; % Ensure it remains empty if no real data
end

%% --- 4. Set up parameters for LFP analysis ---
fprintf('Setting up LFP analysis parameters...\n');

% Validate channel index
if chan < 1 || chan > length(brain_regions) || chan > length(brain_region_names)
     error('Selected channel index "chan" (%d) is out of bounds for brain_regions/brain_region_names lists (size %d).', chan, length(brain_regions));
end
analysis_region_code = brain_regions{chan}; % e.g., '+3-1'
analysis_region_name = brain_region_names{chan}; % e.g., 'STN_L_03_01'

% Channel selection - which LFP channel column to analyze
% Assumes createCombinedTable created columns like TD_key0, TD_key1 etc.
% TD_key<N> corresponds to the (N+1)th channel in timeDomainSettings
td_channel_idx = chan - 1; % 0-based index for TD_key
channel_col_name = sprintf('TD_key%d', td_channel_idx);

if ~ismember(channel_col_name, combinedDataTable.Properties.VariableNames)
    error('Channel column "%s" (derived from chan=%d) not found in combinedDataTable. Available TD columns: %s', ...
          channel_col_name, chan, strjoin(combinedDataTable.Properties.VariableNames(startsWith(combinedDataTable.Properties.VariableNames, 'TD_key')), ', '));
end
fprintf('Analyzing channel: %s (%s)\n', analysis_region_name, channel_col_name);

% Extract time series data
time_vector_ms = combinedDataTable.DerivedTime; % Assumes DerivedTime is POSIX time in MILLISECONDS
signal = combinedDataTable.(channel_col_name);

% --- Get sampling rate with fallback ---
fs = NaN; % Initialize fs as NaN

% Check if timeDomainSettings exists and has a valid samplingRate field
if isfield(loaded_combined, 'timeDomainSettings') && ...
   isstruct(loaded_combined.timeDomainSettings) && ...
   isfield(loaded_combined.timeDomainSettings, 'samplingRate') && ...
   ~isempty(loaded_combined.timeDomainSettings.samplingRate) && ...
   isscalar(loaded_combined.timeDomainSettings.samplingRate) && ... % Ensure it's scalar
   isnumeric(loaded_combined.timeDomainSettings.samplingRate) && ... % Ensure it's numeric
   ~isnan(loaded_combined.timeDomainSettings.samplingRate) && ...    % Ensure it's not NaN
   loaded_combined.timeDomainSettings.samplingRate > 0 % Ensure it's positive

    fs = loaded_combined.timeDomainSettings.samplingRate; % Use the loaded value
    fprintf('Using sampling rate loaded from timeDomainSettings: %.2f Hz\n', fs);

    % Optional sanity check against known rate
    if abs(fs - known_td_sampling_rate_hz) > 1 % Allow small tolerance
         warning('Loaded sampling rate (%.2f Hz) differs significantly from expected rate (%.2f Hz). Check data processing.', fs, known_td_sampling_rate_hz);
         % Decide whether to proceed with loaded or known rate here if desired
         % fs = known_td_sampling_rate_hz; % Uncomment to force known rate if discrepancy found
    end
else
    % If loading failed, issue warning and use hardcoded value
    warning('Could not load valid sampling rate from timeDomainSettings. Using hardcoded value: %d Hz.', known_td_sampling_rate_hz);
    fs = known_td_sampling_rate_hz;

    % Optional: Add the field back to the struct if needed later
    % This might be useful if other parts of the script assume the field exists
    % if isfield(loaded_combined, 'timeDomainSettings') && isstruct(loaded_combined.timeDomainSettings)
    %     loaded_combined.timeDomainSettings(1).samplingRate = fs; % Add it back (use index 1 if it was an array)
    % else
    %     loaded_combined.timeDomainSettings = struct('samplingRate', fs); % Create if missing
    % end
end

% Final check if fs is valid before proceeding
if isnan(fs) || fs <= 0
    error('Failed to determine a valid sampling rate (fs). Cannot proceed.');
end
% --- End sampling rate section ---


% Convert window/step sizes to samples (uses the determined 'fs')
window_size_samples = round(window_size_sec * fs);
step_size_samples_m1 = round(step_size_sec_m1 * fs);

%% --- 5a. METHOD 1: Original method - 60s windows with 5s steps ---
fprintf('\n--- METHOD 1: Calculating aperiodic features (Window: %ds, Step: %ds) ---\n', window_size_sec, step_size_sec_m1);

% Calculate number of windows
num_windows_m1 = floor((length(signal) - window_size_samples) / step_size_samples_m1) + 1;
if num_windows_m1 <= 0
     error('Signal length (%d samples) is too short for the specified window size (%d samples).', length(signal), window_size_samples);
end

% Pre-allocate arrays for results
aperiodic_times_m1_ms = NaN(num_windows_m1, 1); % Store time in ms
aperiodic_exponents_m1 = NaN(num_windows_m1, 1);
aperiodic_offsets_m1 = NaN(num_windows_m1, 1);
aperiodic_r_squared_m1 = NaN(num_windows_m1, 1);
window_med_state_m1 = NaN(num_windows_m1, 1);
window_stim_state_m1 = NaN(num_windows_m1, 1);

% Process each window
fprintf('Processing %d windows for Method 1...\n', num_windows_m1);
for i = 1:num_windows_m1
    % Extract window indices
    start_idx = (i-1) * step_size_samples_m1 + 1;
    end_idx = start_idx + window_size_samples - 1;

    % Ensure indices are within bounds (should be handled by num_windows calculation, but double check)
     if end_idx > length(signal)
         warning('Method 1: Window %d exceeds signal length. Stopping.', i);
         % Trim results if stopped early
         num_windows_m1 = i-1;
         aperiodic_times_m1_ms = aperiodic_times_m1_ms(1:num_windows_m1);
         aperiodic_exponents_m1 = aperiodic_exponents_m1(1:num_windows_m1);
         aperiodic_offsets_m1 = aperiodic_offsets_m1(1:num_windows_m1);
         aperiodic_r_squared_m1 = aperiodic_r_squared_m1(1:num_windows_m1);
         window_med_state_m1 = window_med_state_m1(1:num_windows_m1);
         window_stim_state_m1 = window_stim_state_m1(1:num_windows_m1);
         break;
     end

    % Extract window data and time
    window_data = signal(start_idx:end_idx);
    window_time_ms = time_vector_ms(start_idx:end_idx);

    % Store the center time of this window (in ms)
    aperiodic_times_m1_ms(i) = window_time_ms(round(length(window_time_ms)/2));

    % Calculate window medication and stimulation states (most common state in window)
    % Ensure states vectors are long enough (should be same height as combinedDataTable)
    if end_idx <= length(med_state) && end_idx <= length(stim_state)
        window_med_state_m1(i) = mode(med_state(start_idx:end_idx));
        window_stim_state_m1(i) = mode(stim_state(start_idx:end_idx));
    else
        warning('State vectors shorter than signal index %d. Cannot determine state for window %d.', end_idx, i);
        window_med_state_m1(i) = NaN;
        window_stim_state_m1(i) = NaN;
    end

    % Skip windows with too many NaNs
    nan_frac = sum(isnan(window_data)) / length(window_data);
    if nan_frac > 0.2 % Skip if >20% NaNs
        % NaN values are already pre-allocated
        continue;
    end

    % Handle remaining NaNs (e.g., replace with zero or interpolate)
    % Replacing with zero might slightly affect low-frequency power
    window_data(isnan(window_data)) = 0;

    % Determine appropriate frequency range based on stimulation state
    current_stim_state = window_stim_state_m1(i);
    if ~isnan(current_stim_state) && current_stim_state > 0.5 % State = 1 (ON/High)
        freq_range = stim_freq_range;
    else % State = 0 (OFF/Low) or NaN
        freq_range = med_freq_range;
    end

    % Calculate aperiodic component and store results
    [exponent, offset, r_squared] = calculate_aperiodic_direct(window_data, fs, freq_range, channel_col_name);

    aperiodic_exponents_m1(i) = exponent;
    aperiodic_offsets_m1(i) = offset;
    aperiodic_r_squared_m1(i) = r_squared;

    % Display progress
    if mod(i, 100) == 0 || i == num_windows_m1
        fprintf('  Method 1: Processed %d/%d windows (%.1f%%)\n', i, num_windows_m1, 100*i/num_windows_m1);
    end
end
fprintf('Method 1 processing complete.\n');

% Remove any trailing NaNs if processing stopped early (already handled in loop break)
% Remove windows where FOOOF failed (already pre-filled with NaN)
valid_idx_m1 = ~isnan(aperiodic_exponents_m1);
aperiodic_times_m1_ms = aperiodic_times_m1_ms(valid_idx_m1);
aperiodic_exponents_m1 = aperiodic_exponents_m1(valid_idx_m1);
aperiodic_offsets_m1 = aperiodic_offsets_m1(valid_idx_m1);
aperiodic_r_squared_m1 = aperiodic_r_squared_m1(valid_idx_m1);
window_med_state_m1 = window_med_state_m1(valid_idx_m1);
window_stim_state_m1 = window_stim_state_m1(valid_idx_m1);
fprintf('Method 1: Kept %d valid measurements.\n', length(aperiodic_times_m1_ms));


%% --- 5b. METHOD 2: 60s windows ending exactly at PKG measurement times ---
fprintf('\n--- METHOD 2: Calculating aperiodic features (Window: %ds, Aligned to PKG end times) ---\n', window_size_sec);

% Initialize results arrays
num_pkg_intervals = length(pkg_times); % PKG times are in seconds
aperiodic_times_m2_sec = NaN(num_pkg_intervals, 1); % Store time in seconds
aperiodic_exponents_m2 = NaN(num_pkg_intervals, 1);
aperiodic_offsets_m2 = NaN(num_pkg_intervals, 1);
aperiodic_r_squared_m2 = NaN(num_pkg_intervals, 1);
window_med_state_m2 = NaN(num_pkg_intervals, 1);
window_stim_state_m2 = NaN(num_pkg_intervals, 1);

if isempty(pkg_times)
    fprintf('Method 2: No PKG timestamps available. Skipping calculation.\n');
else
    % Calculate and display the actual step size (typically 120s for PKG data)
    if length(pkg_times) > 1
        pkg_intervals_actual = diff(pkg_times);
        avg_pkg_interval = mean(pkg_intervals_actual);
        fprintf('Method 2: Using %d PKG timestamps with average interval of %.1f seconds.\n', num_pkg_intervals, avg_pkg_interval);
    else
        fprintf('Method 2: Only one PKG timestamp available.\n');
    end

    % Process each window ending at PKG timestamp
    fprintf('Processing %d windows for Method 2...\n', num_pkg_intervals);
    time_vector_sec = time_vector_ms / 1000; % Convert neural time to seconds for comparison

    for i = 1:num_pkg_intervals
        % Calculate window boundaries: window_size_sec window ending exactly at PKG time (in seconds)
        end_time_sec = pkg_times(i);
        start_time_sec = end_time_sec - window_size_sec;

        % Find closest indices in the neural time vector (now in seconds)
        [~, end_idx] = min(abs(time_vector_sec - end_time_sec));
        [~, start_idx] = min(abs(time_vector_sec - start_time_sec));

        % Skip if window is outside the neural data range or invalid
        % Add tolerance (e.g., 5 seconds) for matching times
        if start_idx < 1 || end_idx > length(signal) || start_idx >= end_idx || time_vector_sec(start_idx) > start_time_sec + 5 || time_vector_sec(end_idx) < end_time_sec - 5
            warning('Method 2: Window %d (%.1f - %.1f sec) is outside neural data range or invalid indices [%d %d]. Skipping.', i, start_time_sec, end_time_sec, start_idx, end_idx);
            continue; % Keep NaNs
        end

        % Store the actual end time (PKG time in seconds)
        aperiodic_times_m2_sec(i) = end_time_sec;

        % Extract window data
        window_data = signal(start_idx:end_idx);

        % Calculate window medication and stimulation states (most common state in window)
         if end_idx <= length(med_state) && end_idx <= length(stim_state)
             window_med_state_m2(i) = mode(med_state(start_idx:end_idx));
             window_stim_state_m2(i) = mode(stim_state(start_idx:end_idx));
         else
             warning('State vectors shorter than signal index %d. Cannot determine state for window %d.', end_idx, i);
             window_med_state_m2(i) = NaN;
             window_stim_state_m2(i) = NaN;
         end

        % Skip windows with too many NaNs
        nan_frac = sum(isnan(window_data)) / length(window_data);
        if nan_frac > 0.2 % Skip if >20% NaNs
            continue; % Keep NaNs
        end

        window_data(isnan(window_data)) = 0; % Replace NaNs with zeros

        % Determine appropriate frequency range based on stimulation state
        current_stim_state = window_stim_state_m2(i);
        if ~isnan(current_stim_state) && current_stim_state > 0.5 % State = 1 (ON/High)
            freq_range = stim_freq_range;
        else % State = 0 (OFF/Low) or NaN
            freq_range = med_freq_range;
        end

        % Calculate aperiodic component and store results
        [exponent, offset, r_squared] = calculate_aperiodic_direct(window_data, fs, freq_range, channel_col_name);

        aperiodic_exponents_m2(i) = exponent;
        aperiodic_offsets_m2(i) = offset;
        aperiodic_r_squared_m2(i) = r_squared;

        % Display progress
        if mod(i, 10) == 0 || i == num_pkg_intervals
            fprintf('  Method 2: Processed %d/%d windows (%.1f%%)\n', i, num_pkg_intervals, 100*i/num_pkg_intervals);
        end
    end
    fprintf('Method 2 processing complete.\n');

    % Remove windows where FOOOF failed or window was invalid
    valid_idx_m2 = ~isnan(aperiodic_exponents_m2);
    aperiodic_times_m2_sec = aperiodic_times_m2_sec(valid_idx_m2);
    aperiodic_exponents_m2 = aperiodic_exponents_m2(valid_idx_m2);
    aperiodic_offsets_m2 = aperiodic_offsets_m2(valid_idx_m2);
    aperiodic_r_squared_m2 = aperiodic_r_squared_m2(valid_idx_m2);
    window_med_state_m2 = window_med_state_m2(valid_idx_m2);
    window_stim_state_m2 = window_stim_state_m2(valid_idx_m2);
    fprintf('Method 2: Kept %d valid measurements.\n', length(aperiodic_times_m2_sec));
end


%% --- 5c. METHOD 3: Average Method 1 exponents in 120s window preceding PKG times ---
fprintf('\n--- METHOD 3: Averaging Method 1 results in %ds window preceding PKG times ---\n', pkg_interval_sec);

% Initialize results arrays
aperiodic_times_m3_sec = NaN(num_pkg_intervals, 1); % Store time in seconds
aperiodic_exponents_m3 = NaN(num_pkg_intervals, 1);
aperiodic_offsets_m3 = NaN(num_pkg_intervals, 1);
aperiodic_r_squared_m3 = NaN(num_pkg_intervals, 1);
window_med_state_m3 = NaN(num_pkg_intervals, 1);
window_stim_state_m3 = NaN(num_pkg_intervals, 1);

if isempty(pkg_times)
    fprintf('Method 3: No PKG timestamps available. Skipping calculation.\n');
elseif isempty(aperiodic_times_m1_ms)
     fprintf('Method 3: No valid Method 1 results available to average. Skipping calculation.\n');
else
     % Calculate and display the actual step size (typically 120s for PKG data)
    if length(pkg_times) > 1
        pkg_intervals_actual = diff(pkg_times);
        avg_pkg_interval = mean(pkg_intervals_actual);
        fprintf('Method 3: Using %d PKG timestamps with average interval of %.1f seconds.\n', num_pkg_intervals, avg_pkg_interval);
    else
        fprintf('Method 3: Only one PKG timestamp available.\n');
    end

    % Convert Method 1 times to seconds for comparison
    aperiodic_times_m1_sec = aperiodic_times_m1_ms / 1000;

    % For each PKG timestamp, find all Method 1 windows that fall within the window before it
    fprintf('Processing %d intervals for Method 3...\n', num_pkg_intervals);
    for i = 1:num_pkg_intervals
        % Define the averaging window preceding the PKG timestamp (in seconds)
        pkg_time_sec = pkg_times(i);
        avg_window_start_sec = pkg_time_sec - pkg_interval_sec; % e.g., 120s before PKG time
        avg_window_end_sec = pkg_time_sec;

        % Store the end time (PKG time in seconds)
        aperiodic_times_m3_sec(i) = pkg_time_sec;

        % Find Method 1 results (using their center times in seconds) within this window
        m1_indices_in_window = find(aperiodic_times_m1_sec >= avg_window_start_sec & aperiodic_times_m1_sec < avg_window_end_sec);

        % If no Method 1 windows fall in this interval, skip
        if isempty(m1_indices_in_window)
            continue; % Keep NaNs
        end

        % Average the values within the window
        aperiodic_exponents_m3(i) = mean(aperiodic_exponents_m1(m1_indices_in_window), 'omitnan');
        aperiodic_offsets_m3(i) = mean(aperiodic_offsets_m1(m1_indices_in_window), 'omitnan');
        aperiodic_r_squared_m3(i) = mean(aperiodic_r_squared_m1(m1_indices_in_window), 'omitnan');

        % For states, take the mode (most common state) of the M1 windows within the interval
        window_med_state_m3(i) = mode(window_med_state_m1(m1_indices_in_window));
        window_stim_state_m3(i) = mode(window_stim_state_m1(m1_indices_in_window));

        % Display progress
        if mod(i, 10) == 0 || i == num_pkg_intervals
            fprintf('  Method 3: Processed %d/%d intervals (%.1f%%)\n', i, num_pkg_intervals, 100*i/num_pkg_intervals);
        end
    end
    fprintf('Method 3 processing complete.\n');

    % Remove intervals where averaging failed (e.g., no M1 windows found)
    valid_idx_m3 = ~isnan(aperiodic_exponents_m3);
    aperiodic_times_m3_sec = aperiodic_times_m3_sec(valid_idx_m3);
    aperiodic_exponents_m3 = aperiodic_exponents_m3(valid_idx_m3);
    aperiodic_offsets_m3 = aperiodic_offsets_m3(valid_idx_m3);
    aperiodic_r_squared_m3 = aperiodic_r_squared_m3(valid_idx_m3);
    window_med_state_m3 = window_med_state_m3(valid_idx_m3);
    window_stim_state_m3 = window_stim_state_m3(valid_idx_m3);
    fprintf('Method 3: Kept %d valid measurements.\n', length(aperiodic_times_m3_sec));
end


%% --- 6. Save results from all methods ---
fprintf('\nSaving aperiodic results from all methods...\n');

% Create results tables (Convert POSIX times in seconds/ms back to datetime for saving)
aperiodic_results_m1 = table(datetime(aperiodic_times_m1_ms/1000, 'ConvertFrom', 'posixtime'), ... % Convert ms to sec first
                              aperiodic_exponents_m1, aperiodic_offsets_m1, aperiodic_r_squared_m1, ...
                              window_med_state_m1, window_stim_state_m1, ...
                              'VariableNames', {'Timestamp', 'Exponent', 'Offset', 'R_Squared', 'MedState', 'StimState'});
aperiodic_results_m1.Properties.VariableUnits = {'datetime', '', '', '', '0=OFF/1=ON', '0=LOW/1=HIGH'};
aperiodic_results_m1.Properties.Description = 'Method 1: 60s windows, 5s steps, time is window center';

aperiodic_results_m2 = table(datetime(aperiodic_times_m2_sec, 'ConvertFrom', 'posixtime'), ... % Already in seconds
                              aperiodic_exponents_m2, aperiodic_offsets_m2, aperiodic_r_squared_m2, ...
                              window_med_state_m2, window_stim_state_m2, ...
                              'VariableNames', {'Timestamp', 'Exponent', 'Offset', 'R_Squared', 'MedState', 'StimState'});
aperiodic_results_m2.Properties.VariableUnits = {'datetime', '', '', '', '0=OFF/1=ON', '0=LOW/1=HIGH'};
aperiodic_results_m2.Properties.Description = 'Method 2: 60s windows ending at PKG time';

aperiodic_results_m3 = table(datetime(aperiodic_times_m3_sec, 'ConvertFrom', 'posixtime'), ... % Already in seconds
                              aperiodic_exponents_m3, aperiodic_offsets_m3, aperiodic_r_squared_m3, ...
                              window_med_state_m3, window_stim_state_m3, ...
                              'VariableNames', {'Timestamp', 'Exponent', 'Offset', 'R_Squared', 'MedState', 'StimState'});
aperiodic_results_m3.Properties.VariableUnits = {'datetime', '', '', '', '0=OFF/1=ON', '0=LOW/1=HIGH'};
aperiodic_results_m3.Properties.Description = 'Method 3: Average of Method 1 in 120s window preceding PKG time';


% Create metadata struct
metadata = struct();
metadata.session_id = session_id_to_process;
metadata.device_id = device_id_to_process;
metadata.hemisphere = hemi_to_process;
metadata.analysis_channel_code = analysis_region_code; % e.g., '+3-1'
metadata.analysis_channel_name = analysis_region_name; % e.g., 'STN_L_03_01'
metadata.sampling_rate_hz = fs; % Save the sampling rate actually used
metadata.window_size_sec = window_size_sec;
metadata.med_freq_range = med_freq_range;
metadata.stim_freq_range = stim_freq_range;
metadata.fooof_implementation = 'Direct Python FOOOF via MATLAB-Python interface';
metadata.pkg_interval_sec = pkg_interval_sec;
metadata.pkg_patient_id = pkg_patient_id;
metadata.analysis_date = datetime('now');

% Add method-specific metadata
metadata.m1_description = aperiodic_results_m1.Properties.Description;
metadata.m2_description = aperiodic_results_m2.Properties.Description;
metadata.m3_description = aperiodic_results_m3.Properties.Description;
metadata.m1_step_size_sec = step_size_sec_m1;

% Define base filename for saving
base_output_filename = sprintf('Aperiodic_%s_%s_%s_%s', ...
                               hemi_to_process, session_id_to_process, device_id_to_process, analysis_region_name);

% Save combined file with all methods
results_file_all = fullfile(processed_data_folder, [base_output_filename '_AllMethods.mat']);
save(results_file_all, 'aperiodic_results_m1', 'aperiodic_results_m2', 'aperiodic_results_m3', 'metadata', 'pkg_times', '-v7.3');
fprintf('All methods aperiodic results saved to: %s\n', results_file_all);

% Save individual method files (optional, for specific use cases)
results_file_m1 = fullfile(processed_data_folder, [base_output_filename '_Method1.mat']);
save(results_file_m1, 'aperiodic_results_m1', 'metadata', '-v7.3');
fprintf('Method 1 results saved to: %s\n', results_file_m1);

if ~isempty(pkg_times) % Only save M2/M3 if PKG data was used
    results_file_m2 = fullfile(processed_data_folder, [base_output_filename '_Method2.mat']);
    save(results_file_m2, 'aperiodic_results_m2', 'metadata', '-v7.3');
    fprintf('Method 2 results saved to: %s\n', results_file_m2);

    results_file_m3 = fullfile(processed_data_folder, [base_output_filename '_Method3.mat']);
    save(results_file_m3, 'aperiodic_results_m3', 'metadata', '-v7.3');
    fprintf('Method 3 results saved to: %s\n', results_file_m3);
end

%% --- 7. Visualize and compare the three methods ---
fprintf('Visualizing and comparing methods...\n');

figure('Name', sprintf('Aperiodic Comparison: %s %s %s', hemi_to_process, session_id_to_process, analysis_region_name), ...
       'Position', [100 100 1200 800]);

% Plot aperiodic exponents time series comparison
ax1 = subplot(3,1,1);
hold on; grid on;
if ~isempty(aperiodic_results_m1)
    plot(aperiodic_results_m1.Timestamp, aperiodic_results_m1.Exponent, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
end
if ~isempty(aperiodic_results_m2)
    plot(aperiodic_results_m2.Timestamp, aperiodic_results_m2.Exponent, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
end
if ~isempty(aperiodic_results_m3)
    plot(aperiodic_results_m3.Timestamp, aperiodic_results_m3.Exponent, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg)');
end
title('Aperiodic Exponent Comparison');
ylabel('Exponent');
legend('show', 'Location', 'best');

% Plot aperiodic offset time series comparison
ax2 = subplot(3,1,2);
hold on; grid on;
if ~isempty(aperiodic_results_m1)
    plot(aperiodic_results_m1.Timestamp, aperiodic_results_m1.Offset, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
end
if ~isempty(aperiodic_results_m2)
    plot(aperiodic_results_m2.Timestamp, aperiodic_results_m2.Offset, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
end
if ~isempty(aperiodic_results_m3)
    plot(aperiodic_results_m3.Timestamp, aperiodic_results_m3.Offset, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg)');
end
title('Aperiodic Offset Comparison');
ylabel('Offset');
legend('show', 'Location', 'best');


% Plot goodness of fit comparison
ax3 = subplot(3,1,3);
hold on; grid on;
if ~isempty(aperiodic_results_m1)
    plot(aperiodic_results_m1.Timestamp, aperiodic_results_m1.R_Squared, 'LineWidth', 1, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Method 1 (5s steps)');
end
if ~isempty(aperiodic_results_m2)
    plot(aperiodic_results_m2.Timestamp, aperiodic_results_m2.R_Squared, 'o-', 'LineWidth', 1.5, 'Color', 'r', 'MarkerSize', 4, 'DisplayName', 'Method 2 (aligned to PKG)');
end
if ~isempty(aperiodic_results_m3)
    plot(aperiodic_results_m3.Timestamp, aperiodic_results_m3.R_Squared, 's-', 'LineWidth', 1.5, 'Color', 'b', 'MarkerSize', 4, 'DisplayName', 'Method 3 (120s avg)');
end
title('Goodness of Fit (R^2) Comparison');
ylabel('R^2');
xlabel('Time');
ylim([0 1]); % R^2 should be between 0 and 1
legend('show', 'Location', 'best');


% Link axes for synchronized zooming
linkaxes([ax1, ax2, ax3], 'x');
try % Use dynamicDateTicks if available
    dynamicDateTicks([ax1, ax2, ax3], 'linked');
catch
    warning('dynamicDateTicks not found or failed. Using standard datetick.');
    datetick(ax3, 'x', 'keeplimits'); % Apply standard datetick to the bottom axis
end

% Save figure
fig_filename = fullfile(processed_data_folder, [base_output_filename '_Comparison']);
saveas(gcf, [fig_filename '.fig']);
saveas(gcf, [fig_filename '.png']);
fprintf('Comparison figure saved to: %s (.fig/.png)\n', fig_filename);

% Calculate and display statistics for each method
fprintf('\nMethod statistics (Mean +/- SD):\n');
if ~isempty(aperiodic_results_m1)
    fprintf('  Method 1: %d measurements, Exponent = %.3f +/- %.3f, Offset = %.2f +/- %.2f, R^2 = %.3f +/- %.3f\n', ...
        height(aperiodic_results_m1), ...
        mean(aperiodic_results_m1.Exponent, 'omitnan'), std(aperiodic_results_m1.Exponent, 'omitnan'), ...
        mean(aperiodic_results_m1.Offset, 'omitnan'), std(aperiodic_results_m1.Offset, 'omitnan'), ...
        mean(aperiodic_results_m1.R_Squared, 'omitnan'), std(aperiodic_results_m1.R_Squared, 'omitnan'));
else
     fprintf('  Method 1: No valid measurements.\n');
end
if ~isempty(aperiodic_results_m2)
    fprintf('  Method 2: %d measurements, Exponent = %.3f +/- %.3f, Offset = %.2f +/- %.2f, R^2 = %.3f +/- %.3f\n', ...
         height(aperiodic_results_m2), ...
        mean(aperiodic_results_m2.Exponent, 'omitnan'), std(aperiodic_results_m2.Exponent, 'omitnan'), ...
        mean(aperiodic_results_m2.Offset, 'omitnan'), std(aperiodic_results_m2.Offset, 'omitnan'), ...
        mean(aperiodic_results_m2.R_Squared, 'omitnan'), std(aperiodic_results_m2.R_Squared, 'omitnan'));
else
     fprintf('  Method 2: No valid measurements (or no PKG data).\n');
end
if ~isempty(aperiodic_results_m3)
    fprintf('  Method 3: %d measurements, Exponent = %.3f +/- %.3f, Offset = %.2f +/- %.2f, R^2 = %.3f +/- %.3f\n', ...
         height(aperiodic_results_m3), ...
        mean(aperiodic_results_m3.Exponent, 'omitnan'), std(aperiodic_results_m3.Exponent, 'omitnan'), ...
        mean(aperiodic_results_m3.Offset, 'omitnan'), std(aperiodic_results_m3.Offset, 'omitnan'), ...
        mean(aperiodic_results_m3.R_Squared, 'omitnan'), std(aperiodic_results_m3.R_Squared, 'omitnan'));
else
    fprintf('  Method 3: No valid measurements (or no PKG data / M1 results).\n');
end

fprintf('\nStep 4 processing complete.\n');


%% --- Helper function to calculate aperiodic component using direct Python FOOOF ---
function [exponent, offset, r_squared] = calculate_aperiodic_direct(window_data, fs, freq_range, ~)
    % Default values in case of error
    exponent = NaN;
    offset = NaN;
    r_squared = NaN;

    % Ensure window_data is valid
    if isempty(window_data) || all(window_data == 0) || length(unique(window_data)) < 2
        warning('Window data is empty, all zeros, or constant. Skipping FOOOF calculation.');
        return;
    end
    
    % Handle NaNs
    window_data(isnan(window_data)) = 0;
    
    try
        % Calculate power spectrum using Welch's method
        [pxx, f] = pwelch(window_data, hanning(min(length(window_data), 1024)), [], [], fs);
        
        % Filter to desired frequency range
        freq_mask = f >= freq_range(1) & f <= freq_range(2);
        f_fit = f(freq_mask);
        pxx_fit = pxx(freq_mask);
        
        if length(f_fit) < 6
            warning('Too few frequency points for FOOOF analysis.');
            return;
        end
        
        % Convert to Python arrays (numpy)
        py_freqs = py.numpy.array(f_fit);
        py_spectrum = py.numpy.array(pxx_fit);
        
        % Create FOOOF model with desired settings
        fm = py.fooof.FOOOF(pyargs(...
            'peak_width_limits', py.list({1, 8}), ...
            'max_n_peaks', int32(8), ...
            'min_peak_height', 0.05, ...
            'peak_threshold', 2.0, ...
            'aperiodic_mode', 'fixed'));
        
        % Fit the model
        fm.fit(py_freqs, py_spectrum);
        
        % Extract parameters
        aperiodic_params = double(py.array.array('d', fm.aperiodic_params_));
        if length(aperiodic_params) >= 2
            offset = aperiodic_params(1);
            exponent = aperiodic_params(2);
        end
        
        % Get R-squared value
        r_squared = double(fm.r_squared_);
        
    catch ME
        warning('Error in direct FOOOF calculation: %s', ME.message);
    end
end