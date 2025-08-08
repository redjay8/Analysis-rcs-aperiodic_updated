%% Step2_github.m
%
% This script preprocesses RCS data from a MASTER combinedDataTable (output from Step1_MasterData),
% aligns the CENTER of 120s neural windows with interpolated PKG data timestamps (at 30s intervals),
% calculates PSDs, adds PKG interpolation visualization, and saves the output table to a CSV file
% and metadata to a JSON file.
%
% Assumes Step 1 has generated a 'Step1_MasterData_[Patient]_[Hemi]_AllSessions.mat' file.

clear;
close all;
clc;

%% 1. Configuration Settings
fprintf('Starting Step 2');

script_run_folder = pwd;
%project_base_path = fullfile(script_run_folder, '..'); % Assumes script is run from a 'scripts' or similar folder
project_base_path = '/home/jackson';
% --- Input Step 1  Master Data ---
% <<< ADJUST THE FOLLOWING PARAMETERS AS NEEDED >>>
target_patient_folder_step1 = 'RCS02R'; % fPatient identifier used in Step 1 master filename
target_hemisphere_step1 = 'Right';   % Hemisphere identifier used in Step 1 master filename
% ---
step1_master_mat_filename = sprintf('Step1_MasterData_RCS02R_Right_AllSessions.mat', target_patient_folder_step1, target_hemisphere_step1);
processed_data_folder_step1 = fullfile(project_base_path, 'step1_processed_data_multi_session_final_newnaming_tester'); % Folder where Step 1 master files are saved
load_file_path_step1_master = fullfile(processed_data_folder_step1, step1_master_mat_filename);

if ~isfile(load_file_path_step1_master)
    error('Specified Step 1 Master MAT file not found: %s', load_file_path_step1_master);
end


% --- Toolbox Path ---
toolbox_path = '/home/jackson/Analysis-rcs-data'
fprintf('Analysis-rcs-data toolbox found at: %s\n', toolbox_path);

% --- Output Folder ---
preprocessed_output_folder = fullfile(project_base_path, sprintf('step2_preprocessed_data_120s_neural_aligned_%s_%s_AllSessions_newnaming_tester', target_patient_folder_step1, target_hemisphere_step1));
if ~isfolder(preprocessed_output_folder), mkdir(preprocessed_output_folder); end
fprintf('Step 2 output will be saved to: %s\n', preprocessed_output_folder);

% --- Neural Data Processing Parameters ---
%channels_to_process = {'TD_key0', 'TD_key1', 'TD_key2', 'TD_key3'}; % Channels from combinedDataTable
% stim_blanking_window_ms = 10; % Duration around stim pulse to NaN out
use_stim_blanking = false; % IMPORTANT: Stim blanking on a master combined table is complex.
                         % Assuming false for now, or that it was handled prior to Step 1 combination.
                         % If true, stimLogSettings would need to be aggregated and applied carefully.
% --- Automatically detect all channels of the form Contact_x_x ---

% High-pass filter params - keeping this one since it's actually useful
high_pass_cutoff = 1; % Hz, cutoff frequency for high-pass filter
high_pass_order = 4; % Order for Butterworth high-pass filter design
apply_high_pass_filter = true; % Apply high-pass filter?
min_continuous_chunk_for_filter_sec = 5.0; % Min duration of data without NaNs needed to apply filter

% --- Segmentation and PSD Parameters ---
target_neural_segment_duration_sec = 30.0; % Desired duration of each neural segment for PSD
pwelch_config.window_duration_sec = 2.0; % Increased for better low-freq resolution, e.g., 2.0s or 4.0s
pwelch_config.overlap_percent = 50;
pwelch_config.nfft_multiplier = 1; % NFFT = window length * nfft_multiplier
freq_ranges_for_fooof.LowFreq = [10, 40];
freq_ranges_for_fooof.MidFreq = [30, 90];
freq_ranges_for_fooof.WideFreq = [10, 90];

% --- PKG Data Parameters ---
pkg_interpolation_interval_sec = 30.0;
pkg_base_data_folder = fullfile('/home/jackson', 'PKG'); % <<< ADJUST PATH AS NEEDED
default_pkg_filename_pattern = 'scores_*.csv';

%% 2. Load Step1 Master Processed RC+S Data
fprintf('\nLoading Step1 Master data from: %s\n', load_file_path_step1_master);
data_step1_master = load(load_file_path_step1_master);

required_master_vars = {'combinedDataTable_AllSessions', 'electrode_info_Master', 'all_metaData_Master', 'target_patient_folder', 'target_hemisphere'};
missing_vars = setdiff(required_master_vars, fieldnames(data_step1_master));
if ~isempty(missing_vars)
    error('Missing required variable(s) "%s" from Step1 Master MAT file: %s.', strjoin(missing_vars, ', '), load_file_path_step1_master);
end

combinedDataTable = data_step1_master.combinedDataTable_AllSessions;
electrode_info = data_step1_master.electrode_info_Master;
all_metaData = data_step1_master.all_metaData_Master;


all_variable_names = combinedDataTable.Properties.VariableNames;
contact_channel_pattern = '^key\d+_contact_\d+_\d+$';
channels_to_process = all_variable_names(~cellfun(@isempty, regexp(all_variable_names, contact_channel_pattern, 'once')));
fprintf('Detected %d contact channels to process: %s\n', length(channels_to_process), strjoin(channels_to_process, ', '));


patient_id_from_step1 = data_step1_master.target_patient_folder;
neural_hemisphere = data_step1_master.target_hemisphere;

fprintf('Successfully loaded Master Step1 data for Patient: %s, Hemisphere: %s.\n', patient_id_from_step1, neural_hemisphere);
fprintf('Master combinedDataTable size: %d rows x %d columns.\n', size(combinedDataTable,1), size(combinedDataTable,2));

% Let's figure out the UTC offset situation
if ~isempty(all_metaData) && iscell(all_metaData) && isfield(all_metaData{1}, 'UTCoffset')
    UTCoffset_hours = all_metaData{1}.UTCoffset;
    fprintf('Using UTCoffset = %.1f hours (from the first session in all_metaData_Master).\n', UTCoffset_hours);
    for i_md = 2:length(all_metaData)
        if ~isfield(all_metaData{i_md}, 'UTCoffset') || all_metaData{i_md}.UTCoffset ~= UTCoffset_hours
            warning('UTCoffset varies across sessions in all_metaData_Master. Using %.1f from the first session. This might affect PKG alignment if PKG data spans regions with different offsets represented in the neural data.', UTCoffset_hours);
            break;
        end
    end
elseif ~isempty(all_metaData) && isstruct(all_metaData) && numel(all_metaData) >=1 && isfield(all_metaData(1), 'UTCoffset')
     UTCoffset_hours = all_metaData(1).UTCoffset;
     fprintf('Using UTCoffset = %.1f hours (from the first session in all_metaData_Master struct array).\n', UTCoffset_hours);
else
    warning('UTCoffset not found in the first element of all_metaData_Master or all_metaData_Master is empty. Assuming UTC (offset=0). PKG alignment might be incorrect if local time differs.');
    UTCoffset_hours = 0;
end

%% 3. Load, Prepare, and Interpolate PKG Data (Contralateral)
fprintf('\nLoading and preparing PKG data...\n');

% Figure out which hemisphere we need for PKG (it's the opposite one from neural)
if strcmpi(neural_hemisphere, 'Left')
    contralateral_pkg_hemisphere = 'Right';
elseif strcmpi(neural_hemisphere, 'Right')
    contralateral_pkg_hemisphere = 'Left';
else
    error('Neural hemisphere "%s" is invalid. Must be "Left" or "Right".', neural_hemisphere);
end
fprintf('Neural data is %s hemisphere, looking for PKG data from %s hemisphere.\n', neural_hemisphere, contralateral_pkg_hemisphere);

pkg_hemisphere_folder = fullfile(pkg_base_data_folder, contralateral_pkg_hemisphere);
if ~isfolder(pkg_hemisphere_folder)
    error('PKG data folder for %s hemisphere not found: %s', contralateral_pkg_hemisphere, pkg_hemisphere_folder);
end

pkg_files = dir(fullfile(pkg_hemisphere_folder, default_pkg_filename_pattern));
if isempty(pkg_files)
    error('No PKG score files matching pattern "%s" found in %s.', default_pkg_filename_pattern, pkg_hemisphere_folder);
end
selected_pkg_file = pkg_files(1).name;
fprintf('Found %d PKG file(s) for %s hemisphere, using: %s\n', length(pkg_files), contralateral_pkg_hemisphere, selected_pkg_file);

pkg_file_path = fullfile(pkg_hemisphere_folder, selected_pkg_file);
if ~isfile(pkg_file_path)
    error('Selected PKG file not found: %s', pkg_file_path);
end
fprintf('Loading PKG scores from: %s\n', pkg_file_path);

% Time to wrestle with CSV import options...
opts = detectImportOptions(pkg_file_path);
opts.VariableNamingRule = 'preserve';
opts = setvartype(opts, 'Date_Time', 'string');

pkg_score_cols_to_interp = {'BK', 'DK', 'Tremor_Score', 'Tremor'};
for col_idx = 1:length(pkg_score_cols_to_interp)
    col_name = pkg_score_cols_to_interp{col_idx};
    if ismember(col_name, opts.VariableTypes)
        opts = setvartype(opts, col_name, 'double');
        opts.MissingRule = 'fill';
        opts.FillValue = NaN;
        opts.TreatAsMissing = {'','NA','NaN','N/A','.'};
    end
end
current_warning_state_TAsM = warning('off', 'MATLAB:table:ModifiedAndSpecifiedMissingValues');
pkg_table_orig = readtable(pkg_file_path, opts);
warning(current_warning_state_TAsM);

% Convert any string/cell columns to numeric - gotta make sure everything's numeric
for col_idx = 1:length(pkg_score_cols_to_interp)
    col_name = pkg_score_cols_to_interp{col_idx};
    if ismember(col_name, pkg_table_orig.Properties.VariableNames)
        if iscell(pkg_table_orig.(col_name)) || isstring(pkg_table_orig.(col_name))
            pkg_table_orig.(col_name) = cellfun(@str2doubleq, cellstr(pkg_table_orig.(col_name)));
        elseif ~isnumeric(pkg_table_orig.(col_name))
            try
                pkg_table_orig.(col_name) = double(pkg_table_orig.(col_name));
            catch ME_convert
                 warning('Could not convert PKG column "%s" to double. Setting to NaN. Error: %s', col_name, ME_convert.message);
                 pkg_table_orig.(col_name) = nan(height(pkg_table_orig),1);
            end
        end
    end
end

% BK score needs to be inverted (negative becomes positive, positive becomes zero)
if ismember('BK', pkg_table_orig.Properties.VariableNames) && isnumeric(pkg_table_orig.BK)
    fprintf('Inverting the "BK" score column using element-wise logic.\n');
    positive_indices = pkg_table_orig.BK > 0;
    pkg_table_orig.BK(positive_indices) = 0;
    negative_indices = pkg_table_orig.BK < 0;
    pkg_table_orig.BK(negative_indices) = -pkg_table_orig.BK(negative_indices);
end

if ~ismember('Date_Time', pkg_table_orig.Properties.VariableNames)
    error("PKG table is missing the required 'Date_Time' column.");
end
try
    % Convert PKG times to Unix timestamps
    pkg_datetime_local_naive = datetime(pkg_table_orig.Date_Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', '');
    pkg_datetime_local_aware = pkg_datetime_local_naive;
    pkg_datetime_local_aware.TimeZone = sprintf('%+03d:00', UTCoffset_hours);
    pkg_datetime_utc = datetime(pkg_datetime_local_aware, 'TimeZone', 'UTC');
    pkg_table_orig.PKG_UnixTimestamp_orig = posixtime(pkg_datetime_utc);
catch ME_time
    error('Error converting PKG Date_Time to Unix time. Check format and UTC offset consistency. Error: %s', ME_time.message);
end

% Filter out any "Off_Wrist" data points
if ismember('Off_Wrist', pkg_table_orig.Properties.VariableNames)
    original_height = height(pkg_table_orig);
    if islogical(pkg_table_orig.Off_Wrist)
        pkg_table_orig = pkg_table_orig(~pkg_table_orig.Off_Wrist, :);
    elseif isnumeric(pkg_table_orig.Off_Wrist)
        pkg_table_orig = pkg_table_orig(pkg_table_orig.Off_Wrist == 0, :);
    else
        warning('PKG "Off_Wrist" column has unexpected type. Cannot filter. Check values.');
    end
    fprintf('Filtered PKG data for Off_Wrist. Entries removed: %d. Remaining: %d\n', original_height - height(pkg_table_orig), height(pkg_table_orig));
end

% Prepare table for interpolation
valid_pkg_score_cols = pkg_score_cols_to_interp(ismember(pkg_score_cols_to_interp, pkg_table_orig.Properties.VariableNames));
pkg_table_for_interp = pkg_table_orig(:, [{'PKG_UnixTimestamp_orig', 'Date_Time'}, valid_pkg_score_cols]);
pkg_table_for_interp = sortrows(pkg_table_for_interp, 'PKG_UnixTimestamp_orig');
pkg_table_for_interp = unique(pkg_table_for_interp, 'rows', 'stable');

% Now let's interpolate the PKG data to regular intervals
pkg_table_interp_aligned_scores = table();
if height(pkg_table_for_interp) < 2
    warning('Less than 2 unique, valid PKG data points remain after filtering. Cannot interpolate PKG scores.');
else
    min_pkg_time = min(pkg_table_for_interp.PKG_UnixTimestamp_orig(~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig)));
    max_pkg_time = max(pkg_table_for_interp.PKG_UnixTimestamp_orig(~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig)));

    if isempty(min_pkg_time) || isempty(max_pkg_time)
        warning('No valid PKG timestamps found after filtering. Cannot determine PKG interpolation range.');
    else
        fprintf('Creating target PKG interpolation timestamps every %.1f seconds...\n', pkg_interpolation_interval_sec);
        interp_start_time = ceil(min_pkg_time / pkg_interpolation_interval_sec) * pkg_interpolation_interval_sec;
        interp_end_time = floor(max_pkg_time / pkg_interpolation_interval_sec) * pkg_interpolation_interval_sec;

        if interp_start_time > interp_end_time
             warning('PKG interpolation start time is after end time. Check PKG time range and interval. No PKG interpolation will be performed.');
        else
            PKG_UnixTime_target_interp = (interp_start_time : pkg_interpolation_interval_sec : interp_end_time)';
            pkg_table_interp_aligned_scores = table(PKG_UnixTime_target_interp, 'VariableNames', {'Aligned_PKG_UnixTimestamp'});

            fprintf('Interpolating PKG scores:\n');
            for i_score = 1:length(valid_pkg_score_cols)
                score_name = valid_pkg_score_cols{i_score};
                current_scores_numeric = double(pkg_table_for_interp.(score_name));
                valid_idx_score = ~isnan(current_scores_numeric) & ~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig);

                if sum(valid_idx_score) >= 2
                    unique_times_for_this_score = pkg_table_for_interp.PKG_UnixTimestamp_orig(valid_idx_score);
                    scores_for_this_score = current_scores_numeric(valid_idx_score);
                    [final_unique_times, ia, ic] = unique(unique_times_for_this_score, 'stable');
                    if length(final_unique_times) < length(unique_times_for_this_score)
                        % Average duplicate timepoints
                        final_unique_scores = accumarray(ic, scores_for_this_score, [], @mean);
                    else
                        final_unique_scores = scores_for_this_score(ia);
                    end
                    
                    if length(final_unique_times) >= 2
                        pkg_table_interp_aligned_scores.(['Aligned_' score_name]) = interp1(final_unique_times, final_unique_scores, PKG_UnixTime_target_interp, 'linear', nan);
                        fprintf('  - %s interpolated successfully.\n', score_name);
                    else
                        pkg_table_interp_aligned_scores.(['Aligned_' score_name]) = nan(size(PKG_UnixTime_target_interp));
                        fprintf('  - Warning: Not enough unique non-NaN points (%d) for %s after unique time processing to interpolate.\n', length(final_unique_times), score_name);
                    end
                else
                    pkg_table_interp_aligned_scores.(['Aligned_' score_name]) = nan(size(PKG_UnixTime_target_interp));
                    fprintf('  - Warning: Not enough valid points (sum(valid_idx_score)=%d) to interpolate PKG score: %s\n', sum(valid_idx_score), score_name);
                end
            end
            % Add datetime strings for reference
            if ~isempty(pkg_table_orig) && ~isempty(PKG_UnixTime_target_interp) && ismember('Date_Time', pkg_table_for_interp.Properties.VariableNames) && height(pkg_table_interp_aligned_scores)>0
                [~, nearest_indices] = min(abs(pkg_table_for_interp.PKG_UnixTimestamp_orig' - PKG_UnixTime_target_interp), [], 2);
                pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str = pkg_table_for_interp.Date_Time(nearest_indices);
            elseif height(pkg_table_interp_aligned_scores)>0
                pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str = repmat({''}, height(pkg_table_interp_aligned_scores), 1);
            else
                 pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str = cell(0,1);
            end
        end
    end
end
fprintf('Interpolated PKG data to %d points at %.1fs intervals.\n', height(pkg_table_interp_aligned_scores), pkg_interpolation_interval_sec);
if isempty(pkg_table_interp_aligned_scores) || height(pkg_table_interp_aligned_scores)==0, warning('Interpolated PKG data is empty or has no rows. No segments will be generated.'); end


% %% 3.1 Visualize PKG Interpolation
% fprintf('\nVisualizing PKG interpolation...\n');
% if ~isempty(pkg_table_for_interp) && ~isempty(pkg_table_interp_aligned_scores) && height(pkg_table_for_interp) >= 2 && height(pkg_table_interp_aligned_scores) >=1
%     vis_pkg_output_folder = fullfile(preprocessed_output_folder, 'visualizations_pkg_interp');
%     if ~isfolder(vis_pkg_output_folder), mkdir(vis_pkg_output_folder); end
%     fprintf('Saving PKG interpolation visualizations to: %s\n', vis_pkg_output_folder);
% 
%     for i_score = 1:length(valid_pkg_score_cols)
%         score_name = valid_pkg_score_cols{i_score};
%         aligned_score_name = ['Aligned_' score_name];
% 
%         if ismember(score_name, pkg_table_for_interp.Properties.VariableNames) && ...
%            ismember(aligned_score_name, pkg_table_interp_aligned_scores.Properties.VariableNames)
% 
%             fig_h = figure('Name', ['PKG Interpolation: ' score_name], 'Position', [100, 100, 900, 600], 'Visible', 'off');
%             ax = axes(fig_h);
%             hold(ax, 'on');
% 
%             % Plot original PKG points
%             original_scores_numeric = double(pkg_table_for_interp.(score_name));
%             valid_orig_idx = ~isnan(original_scores_numeric) & ~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig);
%             if sum(valid_orig_idx) > 0
%                 plot(ax, pkg_table_for_interp.PKG_UnixTimestamp_orig(valid_orig_idx), original_scores_numeric(valid_orig_idx), 'o-', 'DisplayName', 'Original PKG points', 'MarkerFaceColor', 'b');
%             end
% 
%             % Plot interpolated points
%             interpolated_scores_numeric = double(pkg_table_interp_aligned_scores.(aligned_score_name));
%             valid_interp_idx = ~isnan(interpolated_scores_numeric) & ~isnan(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
%             if sum(valid_interp_idx) > 0
%                 plot(ax, pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp(valid_interp_idx), interpolated_scores_numeric(valid_interp_idx), '.-', 'DisplayName', ['Interpolated PKG (' num2str(pkg_interpolation_interval_sec) 's)'], 'Color', 'r', 'MarkerSize', 10);
%             end
% 
%             hold(ax, 'off');
%             xlabel(ax, 'Unix Timestamp (seconds)');
%             ylabel(ax, strrep(score_name, '_', ' '));
%             title(ax, ['PKG Interpolation Check: ' strrep(score_name, '_', ' ') ' for ' patient_id_from_step1]);
%             legend(ax, 'show', 'Location', 'best');
%             grid(ax, 'on');
% 
%             % Pretty up the x-axis with actual dates
%             if sum(valid_orig_idx)>0 || sum(valid_interp_idx) >0
%                 xlim_val = xlim(ax);
%                 xticks_val = xticks(ax);
%                 xticklabels_val = cellstr(datetime(xticks_val, 'ConvertFrom', 'posixtime', 'TimeZone',  sprintf('%+03d:00', UTCoffset_hours), 'Format', 'MM-dd HH:mm'));
%                 set(ax, 'XTickLabel', xticklabels_val, 'XTickLabelRotation', 30);
%                 xlabel(ax, sprintf('Time (Local, UTC%+d)', UTCoffset_hours));
%             end
% 
%             try
%                 saveas(fig_h, fullfile(vis_pkg_output_folder, sprintf('PKG_Interp_%s_%s_%s.png', patient_id_from_step1, neural_hemisphere, score_name)));
%                 fprintf('  - Saved visualization for %s\n', score_name);
%             catch ME_save
%                  warning('Could not save PKG visualization for %s. Error: %s', score_name, ME_save.message);
%             end
%             close(fig_h);
%         else
%             fprintf('  - Skipping visualization for %s (column not found in original or interpolated table).\n', score_name);
%         end
%         fprintf('📊 Interpolated PKG points available: %d\n', height(pkg_table_interp_aligned_scores));
%     end
% else
%     fprintf('  Skipping PKG interpolation visualization (insufficient data).\n');
% end
% 

%% 4. Convert Neural Timestamps in combinedDataTable (Master Table)
fprintf('\nVerifying neural data timestamps in master combinedDataTable...\n');
if ~ismember('DerivedTime', combinedDataTable.Properties.VariableNames)
    error('Master combinedDataTable is missing the "DerivedTime" (Unix ms) column.');
end
if ~ismember('localTime', combinedDataTable.Properties.VariableNames) || ~isdatetime(combinedDataTable.localTime)
    error('Master combinedDataTable is missing "localTime" column or it is not a datetime object.');
end

if ~ismember('Neural_UnixTimestamp', combinedDataTable.Properties.VariableNames)
    combinedDataTable.Neural_UnixTimestamp = combinedDataTable.DerivedTime / 1000;
    fprintf('  "Neural_UnixTimestamp" (seconds) column created from "DerivedTime".\n');
else
    fprintf('  "Neural_UnixTimestamp" column already exists. Verifying it is in seconds.\n');
    if abs(median(combinedDataTable.Neural_UnixTimestamp) - median(combinedDataTable.DerivedTime/1000)) > 1
        warning('"Neural_UnixTimestamp" exists but might not be in seconds. Check consistency with DerivedTime (ms).');
    end
end
if ~isnumeric(combinedDataTable.Neural_UnixTimestamp)
    error('"Neural_UnixTimestamp" column is not numeric.');
end

fprintf('\n🧠 Neural & PKG time range check:\n');
min_neural_time = min(combinedDataTable.Neural_UnixTimestamp);
max_neural_time = max(combinedDataTable.Neural_UnixTimestamp);
fprintf('  Neural timestamp range: %.2f to %.2f (seconds)\n', min_neural_time, max_neural_time);

min_pkg_time = min(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
max_pkg_time = max(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
fprintf('  Interpolated PKG timestamp range: %.2f to %.2f (seconds)\n', min_pkg_time, max_pkg_time);


%% 5. Preprocess Neural Channels, Segment, Calculate PSD, Align (using Master Table)
fprintf('\nProcessing neural channels from Master combinedDataTable (using parallel processing)...\n');
fs = 250; % Sampling rate
pwelch_config.fs = fs; % Set fs in pwelch_config here, it will be broadcast to parfor workers
fprintf('Assuming consistent sampling rate Fs = %.1f Hz for all processed channels.\n', fs);

% *** ADD THIS LINE BEFORE PARFOR ***
half_neural_segment_duration_sec = target_neural_segment_duration_sec / 2;
% --- PARFOR MODIFICATION START ---
% Time to unleash the parallel processing power!
num_channels_to_process = length(channels_to_process);
channel_processing_results = cell(1, num_channels_to_process); % Cell array to store results from each worker

% Define this variable BEFORE parfor so it's available for debugging


% %% DEBUG: Check data before parfor processing
% fprintf('\n=== DEBUGGING BEFORE PARFOR ===\n');
% fprintf('PKG interpolated data summary:\n');
% fprintf('  - PKG points: %d\n', height(pkg_table_interp_aligned_scores));
% fprintf('  - PKG time range: %.0f to %.0f\n', min(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp), max(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp));
% 
% fprintf('\nNeural data summary:\n');
% fprintf('  - Neural points: %d\n', height(combinedDataTable));
% fprintf('  - Neural time range: %.0f to %.0f\n', min(combinedDataTable.Neural_UnixTimestamp), max(combinedDataTable.Neural_UnixTimestamp));
% 
% fprintf('\nSegment window parameters:\n');
% fprintf('  - Target segment duration: %.1f seconds\n', target_neural_segment_duration_sec);
% fprintf('  - Half window: %.1f seconds\n', target_neural_segment_duration_sec/2);
% 
% % Test the first few PKG timepoints manually
% fprintf('\nTesting first 3 PKG timepoints for potential segments:\n');
% % half_neural_segment_duration_sec = target_neural_segment_duration_sec / 2;
% for test_i = 1:min(3, height(pkg_table_interp_aligned_scores))
%     pkg_time = pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp(test_i);
%     target_start = pkg_time - half_neural_segment_duration_sec;
%     target_end = pkg_time + half_neural_segment_duration_sec;
%     
%     neural_indices = find(combinedDataTable.Neural_UnixTimestamp >= target_start & ...
%                          combinedDataTable.Neural_UnixTimestamp < target_end);
%     
%     fprintf('  PKG point %d (time=%.0f):\n', test_i, pkg_time);
%     fprintf('    Window: %.0f to %.0f\n', target_start, target_end);
%     fprintf('    Neural points in window: %d\n', length(neural_indices));
%     
%     if ~isempty(neural_indices)
%         % Check first channel as example
%         test_channel = channels_to_process{1};
%         test_data = combinedDataTable.(test_channel)(neural_indices);
%         nan_pct = sum(isnan(test_data)) / length(test_data) * 100;
%         fprintf('    %s: %.1f%% NaNs in this window\n', test_channel, nan_pct);
%     end
% end
% 
% fprintf('\nOverall temporal alignment check:\n');
% overlap_start = max(min(combinedDataTable.Neural_UnixTimestamp), min(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp));
% overlap_end = min(max(combinedDataTable.Neural_UnixTimestamp), max(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp));
% overlap_duration_hours = (overlap_end - overlap_start) / 3600;
% fprintf('  Actual overlap: %.2f hours (%.0f to %.0f)\n', overlap_duration_hours, overlap_start, overlap_end);
% 
% if overlap_duration_hours <= 0
%     fprintf('  ❌ NO ACTUAL OVERLAP! This explains why no segments are generated.\n');
% else
%     fprintf('  ✅ Overlap exists, issue is likely in segment filtering logic.\n');
% end
% 
% fprintf('=== END DEBUGGING ===\n\n');

parfor i_ch = 1:num_channels_to_process
    % Get current channel key for this worker
    channel_key_in_cdt = channels_to_process{i_ch};
    
    % Initialize a local cell array to store segment entries for THIS channel by THIS worker
    segments_for_this_channel = {}; 
    
    % Make local copies of variables that are broadcast or specific to this iteration.
    % This can sometimes improve clarity or avoid unintended sharing issues in complex scenarios.
    % MATLAB's parfor is generally good at handling broadcasting of read-only data.
    local_fs = fs; % fs is defined outside and read-only within the loop
    local_pwelch_config = pwelch_config; % pwelch_config is a struct defined outside
    local_combinedDataTable = combinedDataTable; % For reading channel data and timestamps
    local_electrode_info = electrode_info; % For reading channel labels
    local_pkg_table_interp_aligned_scores = pkg_table_interp_aligned_scores; % Read-only for segmentation

    % Determine electrode label for current channel
    current_electrode_label = channel_key_in_cdt;
    % fprintf inside parfor can be messy due to interleaved output. 
    % Consider removing or logging to a file/variable if detailed per-channel logs are needed.
    % For now, this will print, but expect interleaved messages.
    fprintf('Worker processing Channel: %s (Label: %s)...\n', channel_key_in_cdt, current_electrode_label);

    if ~ismember(channel_key_in_cdt, local_combinedDataTable.Properties.VariableNames)
        % fprintf('  Channel %s not found in combinedDataTable. Worker skipping.\n', channel_key_in_cdt);
        channel_processing_results{i_ch} = {}; % Assign empty cell for this iteration's result
        continue; % Next parfor iteration
    end

    channel_raw_neural_data = local_combinedDataTable.(channel_key_in_cdt);
    channel_raw_unix_time_sec = local_combinedDataTable.Neural_UnixTimestamp;

    % if i_ch == 1  % Only do this once
    %     all_nan_mask = all(ismissing(local_combinedDataTable(:, channels_to_process)), 2);
    %     num_removed = sum(all_nan_mask);
    %     if num_removed > 0
    %         fprintf('Removed %d rows where all selected channels are NaN.\n', num_removed);
    %     end
    %     local_combinedDataTable(all_nan_mask, :) = [];
    % end

    if ~isnumeric(channel_raw_neural_data) || ~isnumeric(channel_raw_unix_time_sec)
        %  fprintf('  Channel %s data or its Unix time is not numeric. Worker skipping.\n', channel_key_in_cdt);
         channel_processing_results{i_ch} = {};
         continue;
    end
    if length(channel_raw_neural_data) ~= length(channel_raw_unix_time_sec)
        % fprintf('  Mismatch between data length and time length for channel %s. Worker skipping.\n', channel_key_in_cdt);
        channel_processing_results{i_ch} = {};
        continue;
    end

    data_processed = double(channel_raw_neural_data);
    
    % DEBUG: Add this for the first channel only to avoid spam
    if i_ch == 1
        fprintf('  DEBUG Channel %s:\n', channel_key_in_cdt);
        fprintf('    Neural data points: %d\n', length(channel_raw_neural_data));
        fprintf('    NaN percentage: %.1f%%\n', sum(isnan(channel_raw_neural_data))/length(channel_raw_neural_data)*100);
        fprintf('    PKG points to check: %d\n', height(local_pkg_table_interp_aligned_scores));
        
        % Test first PKG point
        if height(local_pkg_table_interp_aligned_scores) > 0
            test_pkg_time = local_pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp(1);
            test_start = test_pkg_time - half_neural_segment_duration_sec;
            test_end = test_pkg_time + half_neural_segment_duration_sec;
            test_indices = find(channel_raw_unix_time_sec >= test_start & channel_raw_unix_time_sec < test_end);
            fprintf('    First PKG point test: %d neural samples in window\n', length(test_indices));
            
            if ~isempty(test_indices)
                test_segment = data_processed(test_indices);
                test_nan_pct = sum(isnan(test_segment)) / length(test_segment) * 100;
                fprintf('    Test segment: %d samples, %.1f%% NaNs\n', length(test_segment), test_nan_pct);
            end
        end
    end

    % Apply high-pass filter if requested
    if apply_high_pass_filter
        % fprintf('  Worker %d applying high-pass filter to %s...\n', i_ch, channel_key_in_cdt);
        try
            fir_order = 256;  % You can tune this based on fs and data length
            cutoff_hz = high_pass_cutoff;
            d = designfilt('highpassfir', ...
                'FilterOrder', 256, ...
                'CutoffFrequency', 1, ...
                'SampleRate', fs);

            [data_processed, ~] = applyFilterToValidSegments_FIR(data_processed,d, min_continuous_chunk_for_filter_sec);

            disp(size(data_processed));
        catch 
        end
    end

    % Check if we have PKG data to work with
    if isempty(local_pkg_table_interp_aligned_scores) || height(local_pkg_table_interp_aligned_scores) == 0
        % fprintf('  Worker %d skipping segment generation for %s as interpolated PKG data is empty.\n', i_ch, channel_key_in_cdt);
        channel_processing_results{i_ch} = {};
        continue; 
    end

    % fprintf('  Worker %d generating aligned %.1f-second PSD segments for %s...\n', i_ch, target_neural_segment_duration_sec, channel_key_in_cdt);
    num_segments_this_channel_worker = 0;
    % half_neural_segment_duration_sec = target_neural_segment_duration_sec / 2;

    % Loop through each interpolated PKG timepoint
    for i_pkg_pt = 1:height(local_pkg_table_interp_aligned_scores)
        current_aligned_pkg_unixtime = local_pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp(i_pkg_pt);
        if isnan(current_aligned_pkg_unixtime), continue; end

        % Center the neural window around the PKG timepoint
        target_start_unixtime_sec = current_aligned_pkg_unixtime - half_neural_segment_duration_sec;
        target_end_unixtime_sec = current_aligned_pkg_unixtime + half_neural_segment_duration_sec;

        neural_indices_in_window = find(channel_raw_unix_time_sec >= target_start_unixtime_sec & ...
                                          channel_raw_unix_time_sec < target_end_unixtime_sec);
        
        if isempty(neural_indices_in_window), continue; end

        neural_segment_data = data_processed(neural_indices_in_window);
        time_segment_unix_sec = channel_raw_unix_time_sec(neural_indices_in_window);

        actual_start_unixtime_sec = time_segment_unix_sec(1);
        actual_end_unixtime_datapoint_sec = time_segment_unix_sec(end);
        actual_num_samples = length(neural_segment_data);
        
        % Check if segment duration is close enough to what we want
        expected_duration_span_ideal_sec = target_neural_segment_duration_sec;
        actual_duration_span_sec = actual_end_unixtime_datapoint_sec - actual_start_unixtime_sec;
        % duration_tolerance_sec = max(10/local_fs, 0.005 * target_neural_segment_duration_sec);
        duration_tolerance_sec = 1.5;
        if abs(actual_duration_span_sec - expected_duration_span_ideal_sec) > duration_tolerance_sec, continue; end
        
        % Check for too many NaNs
        nan_percentage = sum(isnan(neural_segment_data)) / actual_num_samples * 100;
        max_nan_percentage = 90;
        if nan_percentage > max_nan_percentage
%             fprintf('[%s] ⚠️ Segment rejected due to %.1f%% NaNs\n', channel_key_in_cdt, nan_percentage);
            continue;
        end
        % Make sure we have enough valid data for pwelch
        valid_points_in_segment = neural_segment_data(~isnan(neural_segment_data));
        win_pwelch_samples = round(local_pwelch_config.window_duration_sec * local_fs);
        if length(valid_points_in_segment) < win_pwelch_samples
            fprintf('[%s] ⚠️ Too few valid points for Welch: %d valid < %d required\n', ...
                channel_key_in_cdt, length(valid_points_in_segment), win_pwelch_samples);
            continue;
        end
        if win_pwelch_samples <= 0, continue; end

        % Calculate PSD using pwelch
        noverlap_pwelch_samples = round(win_pwelch_samples * (local_pwelch_config.overlap_percent / 100));
        nfft_pwelch = win_pwelch_samples * local_pwelch_config.nfft_multiplier;
        
        psd_vector = []; freq_vector = []; % Initialize to ensure they exist
        try
            hanning_window = 0.5 - 0.5 * cos(2 * pi * (0:(win_pwelch_samples - 1))' / (win_pwelch_samples - 1));
            psd_vector = [];
            freq_vector = [];
            [psd_vector, freq_vector] = pwelch(valid_points_in_segment, ...
                                               hanning_window, ...
                                               noverlap_pwelch_samples, ...
                                               nfft_pwelch, ...
                                               local_fs);
        catch ME_pwelch
            fprintf('[%s] PSD error: %s\n', channel_key_in_cdt, ME_pwelch.message);
            continue;
        end

        % Build the output entry struct
        entry = struct();
        entry.PatientID = patient_id_from_step1;
        entry.Hemisphere = neural_hemisphere;
        entry.Channel = channel_key_in_cdt;
        entry.ElectrodeLabel = current_electrode_label;
        entry.Neural_Segment_Start_Unixtime = actual_start_unixtime_sec;
        entry.Neural_Segment_End_Unixtime = target_end_unixtime_sec;
        entry.PSD_Data_Str = strjoin(arrayfun(@(x) sprintf('%.8e',x), psd_vector, 'UniformOutput', false), ';');
        entry.Frequency_Vector_Str = strjoin(arrayfun(@(x) sprintf('%.4f',x), freq_vector, 'UniformOutput', false), ';');
        entry.FS = local_fs;

        % Add aligned PKG data
        aligned_pkg_scores_row = local_pkg_table_interp_aligned_scores(i_pkg_pt, :);
        entry.Aligned_PKG_UnixTimestamp = aligned_pkg_scores_row.Aligned_PKG_UnixTimestamp;
        entry.Aligned_PKG_DateTime_Str = aligned_pkg_scores_row.Aligned_PKG_DateTime_Str{1};

        % Add all the PKG scores
        for score_col_idx = 1:length(valid_pkg_score_cols)
            original_score_name = valid_pkg_score_cols{score_col_idx};
            aligned_score_col_name_in_pkg_table = ['Aligned_' original_score_name];
            output_score_col_name = ['Aligned_' original_score_name];
            if ismember(aligned_score_col_name_in_pkg_table, aligned_pkg_scores_row.Properties.VariableNames)
                 entry.(output_score_col_name) = aligned_pkg_scores_row.(aligned_score_col_name_in_pkg_table);
            else
                entry.(output_score_col_name) = NaN;
            end
        end
        segments_for_this_channel{end+1} = entry; % Add to this worker's list for this channel
        num_segments_this_channel_worker = num_segments_this_channel_worker + 1;


    end 
    
    channel_processing_results{i_ch} = segments_for_this_channel; % Store all segments from this worker
    fprintf('Worker for channel %s finished, generated %d segments.\n', channel_key_in_cdt, num_segments_this_channel_worker);
    
end % --- PARFOR MODIFICATION END ---

% Combine results from all workers after the parfor loop
output_segments_list = {};
total_segments_across_all_channels = 0;
for i_res = 1:num_channels_to_process
    if ~isempty(channel_processing_results{i_res}) && iscell(channel_processing_results{i_res})
        % channel_processing_results{i_res} is a cell array like {entry1, entry2, ...}
        % Concatenate these entries into the main list
        output_segments_list = [output_segments_list, channel_processing_results{i_res}{:}];
        total_segments_across_all_channels = total_segments_across_all_channels + length(channel_processing_results{i_res});
    end
end
fprintf('Total %d aligned %.1f-second PSD segments generated across all channels.\n', total_segments_across_all_channels, target_neural_segment_duration_sec);


%% 6. Save Processed Data as CSV and parameters as JSON
output_csv_filename = sprintf('Step2_Aligned120sPSDs_%s_%s_AllSessions.csv', patient_id_from_step1, neural_hemisphere);
output_csv_path = fullfile(preprocessed_output_folder, output_csv_filename);
output_params_filename = sprintf('Step2_params_120sPSDs_%s_%s_AllSessions.json', patient_id_from_step1, neural_hemisphere);
output_params_path = fullfile(preprocessed_output_folder, output_params_filename);

output_data_table = table(); % Initialize empty table
if ~isempty(output_segments_list)
    try
        % Ensure all structs in output_segments_list have the same fields before converting to table
        % This can be an issue if some channels/segments fail and don't populate all fields
        % A more robust conversion might be needed if field consistency isn't guaranteed
        output_data_table = struct2table([output_segments_list{:}], 'AsArray', true);
    catch ME_struct2table
         error('Error converting output segments to table. Check consistency of fields in "entry". Error: %s. Number of segments: %d', ME_struct2table.message, length(output_segments_list));
    end

    fprintf('\nSaving aligned %.1f-second PSD data (from AllSessions)...\n', target_neural_segment_duration_sec);
    if ~isempty(output_data_table)
        fprintf('CSV Header for %s:\n', output_csv_filename);
        disp(strjoin(output_data_table.Properties.VariableNames, ', '));

        try
            writetable(output_data_table, output_csv_path);
            fprintf('Successfully saved data to CSV: %s\n', output_csv_path);
            fprintf('Total %d aligned %.1f-second segments saved across all channels.\n', height(output_data_table), target_neural_segment_duration_sec);
        catch ME_writetable
            error('Failed to write output table to CSV "%s". Error: %s', output_csv_path, ME_writetable.message);
        end
    else
        warning('Output data table is empty after attempting conversion. CSV file will not be saved.');
    end
else
    warning('No data segments generated. CSV file will not be saved.');
end

% Time to save all our processing parameters for posterity
fprintf('\nSaving processing parameters...\n');
params_for_json.patient_id = patient_id_from_step1;
params_for_json.neural_hemisphere = neural_hemisphere;
params_for_json.pkg_data_hemisphere_used = contralateral_pkg_hemisphere;
%params_for_json.channels_processed = channels_to_process;
%params_for_json.electrode_info_used = electrode_info;
params_for_json.target_neural_segment_duration_sec = target_neural_segment_duration_sec;
params_for_json.pkg_interpolation_interval_sec = pkg_interpolation_interval_sec;
params_for_json.pwelch_config_used_in_matlab = pwelch_config; % Original pwelch_config
params_for_json.freq_ranges_defined_for_fooof = freq_ranges_for_fooof;
params_for_json.preprocessing_steps_applied.high_pass_filter = apply_high_pass_filter;
if apply_high_pass_filter, params_for_json.preprocessing_details.high_pass_cutoff_Hz = high_pass_cutoff; end
params_for_json.preprocessing_steps_applied.stim_blanking = use_stim_blanking;
params_for_json.source_step1_master_file = step1_master_mat_filename;
params_for_json.source_pkg_file_used = selected_pkg_file;
params_for_json.utc_offset_hours_used = UTCoffset_hours;
params_for_json.alignment_logic_description = sprintf('Center Alignment: Center of %.1fs neural segment aligned with Aligned_PKG_UnixTimestamp (Window=[T_pkg-%.1fs, T_pkg+%.1fs))', target_neural_segment_duration_sec, target_neural_segment_duration_sec/2, target_neural_segment_duration_sec/2);
params_for_json.script_name = mfilename('fullpath');
params_for_json.script_last_run_timestamp_utc = datetime('now', 'TimeZone', 'UTC', 'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSZ');
params_for_json.matlab_version = version;
params_for_json.step2_output_folder = preprocessed_output_folder;
params_for_json.num_segments_generated = height(output_data_table); % Use height from final table


try
    json_string = jsonencode(params_for_json, 'PrettyPrint', true);
    fid = fopen(output_params_path, 'w');
    if fid == -1
        error('Cannot create JSON parameter file for writing: %s', output_params_path);
    end
    fprintf(fid, '%s', json_string);
    fclose(fid);
    fprintf('Parameters saved to JSON: %s\n', output_params_path);
catch ME_json
    warning('Could not save parameters to JSON file (%s). Error: %s. Saving as MAT instead.', output_params_path, ME_json.message);
    output_params_mat_path = strrep(output_params_path, '.json', '.mat');
    try
        save(output_params_mat_path, 'params_for_json');
        fprintf('Parameters saved to MAT: %s\n', output_params_mat_path);
    catch ME_savemat
         warning('Could not save parameters to MAT file either (%s). Error: %s', output_params_mat_path, ME_savemat.message);
    end
end

fprintf('\nStep 2 (%.1fs Neural Windows Version, using Master Step1 data, with parallel processing) completed.\n', target_neural_segment_duration_sec);

%% ---- Helper Function: applyFilterToValidSegments ----
function [filtered_data, segments_too_short] = applyFilterToValidSegments(data, fs_filt, b, a, min_duration_sec_for_filter) % Renamed fs to fs_filt to avoid conflict if fs is global
    % This function applies a filter to continuous non-NaN segments of data
    % It's smart enough to skip segments that are too short for filtering
    
    if ~isvector(data)
        error('Input data must be a vector.');
    end
    if ~isnumeric(data) || ~isa(data,'double')
       data = double(data);
    end

    % Find where the NaN segments are
    nan_indices = isnan(data);
    transitions = diff([1; nan_indices(:); 1]);
    segment_starts = find(transitions == -1);
    segment_ends = find(transitions == 1) - 1;

    filtered_data = data;
    segments_too_short = 0;

    % Figure out minimum samples needed for filtfilt
    filter_order_for_filtfilt = max(length(b), length(a)) -1;
    min_samples_for_filtfilt = 3 * filter_order_for_filtfilt;
    if min_samples_for_filtfilt < 10, min_samples_for_filtfilt = 10; end

    min_samples_duration_based = round(min_duration_sec_for_filter * fs_filt);
    min_total_samples_required = max(min_samples_for_filtfilt, min_samples_duration_based);
    
    % Handle case where there are no NaN segments (all data is continuous)
    if isempty(segment_starts) && ~all(nan_indices)
        if length(data) < min_total_samples_required
            segments_too_short = 1;
            filtered_data(:) = NaN;
            % warning('Entire dataset is too short for filtering and set to NaN.'); % Suppress in parfor
        else
            try
                filtered_data = filtfilt(b, a, data);
            catch ME_filtfilt_all
                % warning('Filtering error on entire dataset. Data set to NaN. Error: %s', ME_filtfilt_all.message); % Suppress in parfor
                filtered_data(:) = NaN;
            end
        end
        return;
    elseif all(nan_indices)
        % Everything is NaN, nothing to filter
        return;
    end

    % Filter each valid segment
    for i = 1:length(segment_starts)
        seg_start_idx = segment_starts(i);
        seg_end_idx = segment_ends(i);
        current_segment_length = seg_end_idx - seg_start_idx + 1;

        if current_segment_length < min_total_samples_required
            % Segment too short, can't filter it properly
            segments_too_short = segments_too_short + 1;
            filtered_data(seg_start_idx:seg_end_idx) = NaN;
            continue;
        end

        segment_data_for_filt = data(seg_start_idx:seg_end_idx);
        
        try
            segment_filtered = filtfilt(b, a, segment_data_for_filt);
            filtered_data(seg_start_idx:seg_end_idx) = segment_filtered;
        catch ME_filtfilt
            % warning('Filtering error on segment %d. Segment set to NaN. Error: %s', i, ME_filtfilt.message); % Suppress in parfor
            filtered_data(seg_start_idx:seg_end_idx) = NaN;
        end
    end
end

%% ---- Helper Function: str2doubleq ----
function out_val = str2doubleq(str_val_cell_element)
    % Quick and dirty string to double conversion that handles edge cases
    if isempty(str_val_cell_element)
        out_val = NaN;
    elseif isnumeric(str_val_cell_element)
         out_val = double(str_val_cell_element);
    elseif ischar(str_val_cell_element) || isstring(str_val_cell_element)
        str_val = strtrim(char(str_val_cell_element));
        if isempty(str_val)
            out_val = NaN;
        else
            num_val = str2double(str_val);
            if isempty(num_val) && ~isnan(str2double(str_val)) % check if it was non-numeric but not 'NaN' string
                out_val = NaN;
            else
                out_val = num_val;
            end
        end
    else
         out_val = NaN;
    end
end


function [filtered_data, segments_too_short] = applyFilterToValidSegments_FIR(data, d, min_duration_sec_for_filter)
    % Applies FIR filter object d to non-NaN segments using filtfilt

    fs_filt = d.SampleRate;
    if ~isvector(data)
        error('Input data must be a vector.');
    end
    if ~isnumeric(data) || ~isa(data,'double')
        data = double(data);
    end

    nan_indices = isnan(data);
    transitions = diff([1; nan_indices(:); 1]);
    segment_starts = find(transitions == -1);
    segment_ends = find(transitions == 1) - 1;

    filtered_data = data;
    segments_too_short = 0;

    min_samples_duration_based = round(min_duration_sec_for_filter * fs_filt);
    min_samples_for_filtfilt = 3 * d.FilterOrder;  % FilterOrder works for FIR
    min_total_samples_required = max(min_samples_for_filtfilt, min_samples_duration_based);

    if isempty(segment_starts) && ~all(nan_indices)
        if length(data) < min_total_samples_required
            segments_too_short = 1;
            filtered_data(:) = NaN;
        else
            try
                filtered_data = filtfilt(d, data);
            catch
                filtered_data(:) = NaN;
            end
        end
        return;
    elseif all(nan_indices)
        return;
    end

    for i = 1:length(segment_starts)
        seg_start_idx = segment_starts(i);
        seg_end_idx = segment_ends(i);
        seg_len = seg_end_idx - seg_start_idx + 1;

        if seg_len < min_total_samples_required
            segments_too_short = segments_too_short + 1;
            filtered_data(seg_start_idx:seg_end_idx) = NaN;
            continue;
        end

        segment = data(seg_start_idx:seg_end_idx);
        try
            filtered_segment = filtfilt(d, segment);
            filtered_data(seg_start_idx:seg_end_idx) = filtered_segment;
        catch
            filtered_data(seg_start_idx:seg_end_idx) = NaN;
        end
    end
end
