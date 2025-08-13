%% Step2_github.m  —  PARALLELIZED (parfor over PKG epochs)
%
% Preprocesses RCS data from Step 1 Master, aligns CENTER of neural windows
% with interpolated PKG times, computes PSDs, and saves CSV + JSON.
%
% Key parallel change:
%   * parfor iterates over PKG epochs (thousands) instead of channels (few).
%   * Heavy read-only data shipped once per worker via arallel.pool.Constant.

clear; close all; clc;

%% 1. Configuration Settings
fprintf('Starting Step 2\n');

script_run_folder = pwd;
project_base_path = '/home/jackson';

% --- Input Step 1 Master Data ---
target_patient_folder_step1  = 'RCS02R';
target_hemisphere_step1      = 'Right';
step1_master_mat_filename    = sprintf('Step1_MasterData_%s_%s_AllSessions.mat', ...
                                       target_patient_folder_step1, target_hemisphere_step1);
processed_data_folder_step1  = fullfile(project_base_path, 'step1_processed_data_multi_session_final_new');
load_file_path_step1_master  = fullfile(processed_data_folder_step1, step1_master_mat_filename);
if ~isfile(load_file_path_step1_master)
    error('Specified Step 1 Master MAT file not found: %s', load_file_path_step1_master);
end

% --- Toolbox Path ---
toolbox_path = '/home/jackson/Analysis-rcs-data';
fprintf('Analysis-rcs-data toolbox found at: %s\n', toolbox_path);
addpath(genpath(toolbox_path));

% --- Output Folder ---
preprocessed_output_folder = fullfile(project_base_path, ...
    sprintf('step2_preprocessed_data_120s_neural_aligned_%s_%s_AllSessions_newnaming_tester_left_neww', ...
    target_patient_folder_step1, target_hemisphere_step1));
if ~isfolder(preprocessed_output_folder), mkdir(preprocessed_output_folder); end
fprintf('Step 2 output will be saved to: %s\n', preprocessed_output_folder);

% --- PKG Data Parameters ---
pkg_interpolation_interval_sec = 30.0;
pkg_base_data_folder = fullfile('/home/jackson', 'PKG');
default_pkg_filename_pattern = 'scores_*.csv';

% --- Neural Preprocessing ---
apply_high_pass_filter = true;     % your original flag (set as you wish)
high_pass_cutoff       = 1;        % Hz
min_continuous_chunk_for_filter_sec = 5.0;

% --- Segmentation and PSD Parameters ---
target_neural_segment_duration_sec = 30.0;   % your current setting
pwelch_config.window_duration_sec = 2.0;
pwelch_config.overlap_percent     = 50;
pwelch_config.nfft_multiplier     = 1;

% --- FOOOF ranges (metadata only here) ---
freq_ranges_for_fooof.LowFreq  = [10, 40];
freq_ranges_for_fooof.MidFreq  = [30, 90];
freq_ranges_for_fooof.WideFreq = [10, 90];

% --- Stim blanking flag (not used here) ---
use_stim_blanking = false;

%% 2. Load Step1 Master Processed RC+S Data
fprintf('\nLoading Step1 Master data from: %s\n', load_file_path_step1_master);
data_step1_master = load(load_file_path_step1_master);

required_master_vars = {'combinedDataTable_AllSessions','electrode_info_Master','all_metaData_Master','target_patient_folder','target_hemisphere'};
missing_vars = setdiff(required_master_vars, fieldnames(data_step1_master));
if ~isempty(missing_vars)
    error('Missing required variable(s) "%s" from Step1 Master MAT file: %s.', strjoin(missing_vars, ', '), load_file_path_step1_master);
end

combinedDataTable = data_step1_master.combinedDataTable_AllSessions;
electrode_info    = data_step1_master.electrode_info_Master; %#ok<NASGU>  % not needed in loop, kept for metadata
all_metaData      = data_step1_master.all_metaData_Master;
patient_id_from_step1 = data_step1_master.target_patient_folder;
neural_hemisphere     = data_step1_master.target_hemisphere;

% Detect channels like 'keyX_contact_Y_Z'
all_variable_names = combinedDataTable.Properties.VariableNames;
contact_channel_pattern = '^key\d+_contact_\d+_\d+$';
channels_to_process = all_variable_names(~cellfun(@isempty, regexp(all_variable_names, contact_channel_pattern, 'once')));
fprintf('Detected %d contact channels to process: %s\n', length(channels_to_process), strjoin(channels_to_process, ', '));

% UTC offset
if ~isempty(all_metaData) && iscell(all_metaData) && isfield(all_metaData{1}, 'UTCoffset')
    UTCoffset_hours = all_metaData{1}.UTCoffset;
elseif ~isempty(all_metaData) && isstruct(all_metaData) && numel(all_metaData)>=1 && isfield(all_metaData(1),'UTCoffset')
    UTCoffset_hours = all_metaData(1).UTCoffset;
else
    warning('UTCoffset not found; assuming 0.');
    UTCoffset_hours = 0;
end

%% 3. Load, Prepare, Interpolate PKG Data (Contralateral)
fprintf('\nLoading and preparing PKG data...\n');
if strcmpi(neural_hemisphere, 'Left'), contralateral_pkg_hemisphere = 'Right';
elseif strcmpi(neural_hemisphere, 'Right'), contralateral_pkg_hemisphere = 'Left';
else, error('Invalid neural hemisphere: %s', neural_hemisphere);
end
pkg_hemisphere_folder = fullfile(pkg_base_data_folder, contralateral_pkg_hemisphere);
if ~isfolder(pkg_hemisphere_folder), error('PKG folder not found: %s', pkg_hemisphere_folder); end
pkg_files = dir(fullfile(pkg_hemisphere_folder, default_pkg_filename_pattern));
if isempty(pkg_files), error('No PKG files matching %s in %s', default_pkg_filename_pattern, pkg_hemisphere_folder); end
selected_pkg_file = pkg_files(1).name;
pkg_file_path = fullfile(pkg_hemisphere_folder, selected_pkg_file);
fprintf('Loading PKG scores from: %s\n', pkg_file_path);

opts = detectImportOptions(pkg_file_path);
opts.VariableNamingRule = 'preserve';

% Keep Date_Time as string for robust parsing
opts = setvartype(opts, 'Date_Time', 'string');

% Numeric PKG columns we care about
pkg_score_cols_to_interp = {'BK','DK','Tremor_Score','Tremor'};
valid_cols = pkg_score_cols_to_interp(ismember(pkg_score_cols_to_interp, opts.VariableNames));

% Treat common empties as missing BEFORE typing to double (version-safe)
if ~isempty(valid_cols)
    opts = setvaropts(opts, valid_cols, 'TreatAsMissing', {'', 'NA', 'NaN', 'N/A', '.'});
    % If your CSV sometimes uses thousands separators like "1,234"
    try
        opts = setvaropts(opts, valid_cols, 'ThousandsSeparator', ',');
    catch
        % Older MATLAB may not support ThousandsSeparator; safe to ignore
    end
    % Now set numeric type
    for i = 1:numel(valid_cols)
        opts = setvartype(opts, valid_cols{i}, 'double');
    end
end

% Read the table (no FillValue / MissingRule usage)
pkg_table_orig = readtable(pkg_file_path, opts);

% Post-read coercion guard (handles any stragglers)
for i = 1:numel(valid_cols)
    c = valid_cols{i};
    if ~isnumeric(pkg_table_orig.(c))
        pkg_table_orig.(c) = str2double(string(pkg_table_orig.(c)));
    end
end


% Normalise BK (your original logic)
if ismember('BK', pkg_table_orig.Properties.VariableNames) && isnumeric(pkg_table_orig.BK)
    pos = pkg_table_orig.BK > 0; pkg_table_orig.BK(pos) = 0;
    neg = pkg_table_orig.BK < 0; pkg_table_orig.BK(neg) = -pkg_table_orig.BK(neg);
end

if ~ismember('Date_Time', pkg_table_orig.Properties.VariableNames)
    error("PKG table missing 'Date_Time'.");
end
try
    pkg_dt_local_naive = datetime(pkg_table_orig.Date_Time, 'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','');
    pkg_dt_local_aware = pkg_dt_local_naive; pkg_dt_local_aware.TimeZone = sprintf('%+03d:00', UTCoffset_hours);
    pkg_dt_utc = datetime(pkg_dt_local_aware, 'TimeZone','UTC');
    pkg_table_orig.PKG_UnixTimestamp_orig = posixtime(pkg_dt_utc);
catch ME_time
    error('PKG time conversion failed: %s', ME_time.message);
end

% Filter Off_Wrist if present
if ismember('Off_Wrist', pkg_table_orig.Properties.VariableNames)
    H0 = height(pkg_table_orig);
    if islogical(pkg_table_orig.Off_Wrist)
        pkg_table_orig = pkg_table_orig(~pkg_table_orig.Off_Wrist,:);
    elseif isnumeric(pkg_table_orig.Off_Wrist)
        pkg_table_orig = pkg_table_orig(pkg_table_orig.Off_Wrist==0,:);
    end
    fprintf('Filtered Off_Wrist: removed %d rows.\n', H0 - height(pkg_table_orig));
end

% Prepare for interpolation
valid_cols = pkg_score_cols_to_interp(ismember(pkg_score_cols_to_interp, pkg_table_orig.Properties.VariableNames));
pkg_table_for_interp = pkg_table_orig(:, [{'PKG_UnixTimestamp_orig','Date_Time'}, valid_cols]);
pkg_table_for_interp = sortrows(pkg_table_for_interp, 'PKG_UnixTimestamp_orig');
pkg_table_for_interp = unique(pkg_table_for_interp, 'rows', 'stable');

pkg_table_interp_aligned_scores = table();
if height(pkg_table_for_interp) >= 2
    min_pkg_time = min(pkg_table_for_interp.PKG_UnixTimestamp_orig(~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig)));
    max_pkg_time = max(pkg_table_for_interp.PKG_UnixTimestamp_orig(~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig)));
    if ~isempty(min_pkg_time) && ~isempty(max_pkg_time)
        fprintf('Creating PKG timestamps every %.1fs...\n', pkg_interpolation_interval_sec);
        t0 = ceil(min_pkg_time/pkg_interpolation_interval_sec)*pkg_interpolation_interval_sec;
        t1 = floor(max_pkg_time/pkg_interpolation_interval_sec)*pkg_interpolation_interval_sec;
        if t0 <= t1
            PKG_UnixTime_target_interp = (t0:pkg_interpolation_interval_sec:t1)';
            pkg_table_interp_aligned_scores = table(PKG_UnixTime_target_interp, 'VariableNames', {'Aligned_PKG_UnixTimestamp'});
            fprintf('Interpolating PKG scores:\n');
            for i=1:numel(valid_cols)
                nm = valid_cols{i};
                x = double(pkg_table_for_interp.(nm));
                v = ~isnan(x) & ~isnan(pkg_table_for_interp.PKG_UnixTimestamp_orig);
                if nnz(v) >= 2
                    tt = pkg_table_for_interp.PKG_UnixTimestamp_orig(v);
                    yy = x(v);
                    [tt_u, ia, ic] = unique(tt, 'stable');
                    if numel(tt_u) < numel(tt)
                        yy_u = accumarray(ic, yy, [], @mean);
                    else
                        yy_u = yy(ia);
                    end
                    if numel(tt_u) >= 2
                        pkg_table_interp_aligned_scores.(['Aligned_' nm]) = interp1(tt_u, yy_u, PKG_UnixTime_target_interp, 'linear', nan);
                        fprintf('  - %s OK\n', nm);
                    else
                        pkg_table_interp_aligned_scores.(['Aligned_' nm]) = nan(size(PKG_UnixTime_target_interp));
                        fprintf('  - %s skipped (not enough unique pts)\n', nm);
                    end
                else
                    pkg_table_interp_aligned_scores.(['Aligned_' nm]) = nan(size(PKG_UnixTime_target_interp));
                    fprintf('  - %s skipped (not enough valid pts)\n', nm);
                end
            end
            if height(pkg_table_interp_aligned_scores)>0
                [~, idx_near] = min(abs(pkg_table_for_interp.PKG_UnixTimestamp_orig' - PKG_UnixTime_target_interp), [], 2);
                pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str = pkg_table_for_interp.Date_Time(idx_near);
            end
        end
    end
end
fprintf('Interpolated PKG points: %d at %.1fs intervals.\n', height(pkg_table_interp_aligned_scores), pkg_interpolation_interval_sec);

%% 4. Neural timestamps check
fprintf('\nVerifying master combinedDataTable timestamps...\n');
if ~ismember('DerivedTime', combinedDataTable.Properties.VariableNames)
    error('Missing DerivedTime (ms) in master table.');
end
if ~ismember('localTime', combinedDataTable.Properties.VariableNames) || ~isdatetime(combinedDataTable.localTime)
    error('Missing or invalid localTime in master table.');
end
if ~ismember('Neural_UnixTimestamp', combinedDataTable.Properties.VariableNames)
    combinedDataTable.Neural_UnixTimestamp = combinedDataTable.DerivedTime / 1000;
    fprintf('  Created Neural_UnixTimestamp from DerivedTime.\n');
end
if ~isnumeric(combinedDataTable.Neural_UnixTimestamp)
    error('Neural_UnixTimestamp is not numeric.');
end

fprintf('\n🧠 Neural & PKG time range check:\n');
min_neural_time = min(combinedDataTable.Neural_UnixTimestamp);
max_neural_time = max(combinedDataTable.Neural_UnixTimestamp);
fprintf('  Neural: %.2f to %.2f (s)\n', min_neural_time, max_neural_time);
if ~isempty(pkg_table_interp_aligned_scores)
    min_pkg_time = min(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
    max_pkg_time = max(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
    fprintf('  PKG:    %.2f to %.2f (s)\n', min_pkg_time, max_pkg_time);
end

%% 5. PARALLEL PROCESSING (parfor over PKG epochs)
fprintf('\nProcessing neural PSDs with parallelization over PKG epochs...\n');

% === Parallel pool ===
try
    c = parcluster('local');
    if isempty(gcp('nocreate')), parpool(c, c.NumWorkers); end
catch ME
    warning('Could not start parpool: %s. Continuing in serial.', ME.message);
end
% Optional: prevent oversubscription if MKL is multithreaded
% maxNumCompThreads(1);

% Welch params (precompute once)
fs = 250;  pwelch_config.fs = fs;
win_samples = round(pwelch_config.window_duration_sec * fs);
win_samples = max(4, win_samples);
han = 0.5 - 0.5 * cos(2*pi*(0:(win_samples-1))'/(win_samples-1));
noverlap = round(win_samples * (pwelch_config.overlap_percent/100));
nfft     = win_samples * pwelch_config.nfft_multiplier;
welch_const = parallel.pool.Constant(struct('w',han,'noverlap',noverlap,'nfft',nfft,'fs',fs,'win_samples',win_samples));

% High-pass filter (design once)
if apply_high_pass_filter
    d = designfilt('highpassfir','FilterOrder',256,'CutoffFrequency',high_pass_cutoff,'SampleRate',fs);
    filt_const = parallel.pool.Constant(d);
else
    filt_const = [];
end

% Hoist time vector and channel arrays (filter once per channel if enabled)
neural_time_vec   = combinedDataTable.Neural_UnixTimestamp;
neural_time_const = parallel.pool.Constant(neural_time_vec);

num_channels_to_process = numel(channels_to_process);
chan_arrays = cell(1, num_channels_to_process);
for k = 1:num_channels_to_process
    v = combinedDataTable.(channels_to_process{k});
    v = double(v);
    if apply_high_pass_filter
        try
            [v, ~] = applyFilterToValidSegments_FIR(v, filt_const.Value, min_continuous_chunk_for_filter_sec);
        catch
            % keep unfiltered on error
        end
    end
    chan_arrays{k} = v;
end
chan_const = parallel.pool.Constant(chan_arrays);

% Compact PKG struct, restricted to overlapping time window
half_win = target_neural_segment_duration_sec/2;
pkg_struct = struct();
if ~isempty(pkg_table_interp_aligned_scores)
    pkg_t_all = pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp;
    valid_pkg = ~isnan(pkg_t_all) & pkg_t_all >= (min_neural_time+half_win) & pkg_t_all <= (max_neural_time-half_win);
    pkg_struct.t  = pkg_t_all(valid_pkg);
    if ismember('Aligned_PKG_DateTime_Str', pkg_table_interp_aligned_scores.Properties.VariableNames)
        pkg_struct.dt = pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str(valid_pkg);
    else
        pkg_struct.dt = repmat({''}, nnz(valid_pkg), 1);
    end
    for i=1:length(valid_cols)
        nm = ['Aligned_' valid_cols{i}];
        if ismember(nm, pkg_table_interp_aligned_scores.Properties.VariableNames)
            pkg_struct.(nm) = pkg_table_interp_aligned_scores.(nm)(valid_pkg);
        else
            pkg_struct.(nm) = nan(nnz(valid_pkg),1);
        end
    end
else
    pkg_struct.t  = [];
    pkg_struct.dt = {};
end
pkg_const = parallel.pool.Constant(pkg_struct);

% Loop constants
max_nan_pct = 90;
tol_sec     = 1.5;

% === parfor over PKG epochs ===
nPKG = numel(pkg_const.Value.t);
seg_results = cell(nPKG,1);

parfor ip = 1:nPKG
    % Local views
    t   = neural_time_const.Value;
    wc  = welch_const.Value;
    pkg = pkg_const.Value;

    center_t = pkg.t(ip);
    start_t  = center_t - half_win;
    end_t    = center_t + half_win;

    i1 = find(t >= start_t, 1, 'first');
    i2 = find(t <  end_t,  1, 'last');
    if isempty(i1) || isempty(i2) || i2 < i1
        seg_results{ip} = {};
        continue;
    end

    seg_entries = cell(1, num_channels_to_process);
    seg_t = t(i1:i2);
    dur = seg_t(end) - seg_t(1);
    if abs(dur - 2*half_win) > tol_sec
        seg_results{ip} = {};
        continue;
    end

    for ch = 1:num_channels_to_process
        data = chan_const.Value{ch};
        seg  = data(i1:i2);

        nan_pct = 100*sum(isnan(seg))/numel(seg);
        if nan_pct > max_nan_pct, continue; end

        vseg = seg(~isnan(seg));
        if numel(vseg) < wc.win_samples, continue; end

        try
            [Pxx,F] = pwelch(vseg, wc.w, wc.noverlap, wc.nfft, wc.fs);
        catch
            continue;
        end

        e = struct();
        e.PatientID  = patient_id_from_step1;
        e.Hemisphere = neural_hemisphere;
        e.Channel    = channels_to_process{ch};
        e.ElectrodeLabel = e.Channel;
        e.Neural_Segment_Start_Unixtime = seg_t(1);
        e.Neural_Segment_End_Unixtime   = end_t;

        % Keep CSV compatibility (stringify here). If needed, move stringify after parfor.
        e.PSD_Data_Str         = strjoin(compose('%.8e', Pxx), ';');
        e.Frequency_Vector_Str = strjoin(compose('%.4f', F),   ';');
        e.FS = wc.fs;

        e.Aligned_PKG_UnixTimestamp = center_t;
        if ~isempty(pkg.dt), e.Aligned_PKG_DateTime_Str = pkg.dt{ip}; else, e.Aligned_PKG_DateTime_Str = ''; end
        for sc = 1:length(valid_cols)
            nm = ['Aligned_' valid_cols{sc}];
            e.(nm) = pkg.(nm)(ip);
        end
        seg_entries{ch} = e;
    end

    seg_entries = seg_entries(~cellfun('isempty', seg_entries));
    seg_results{ip} = seg_entries;
end

% Collect to one list
output_segments_list = [seg_results{:}];
output_segments_list = output_segments_list(~cellfun('isempty', output_segments_list));
fprintf('Total %d aligned %.1f-second PSD segments generated.\n', numel(output_segments_list), target_neural_segment_duration_sec);

%% 6. Save CSV and JSON
output_csv_filename  = sprintf('Step2_Aligned120sPSDs_%s_%s_AllSessions.csv', patient_id_from_step1, neural_hemisphere);
output_csv_path      = fullfile(preprocessed_output_folder, output_csv_filename);
output_params_filename = sprintf('Step2_params_120sPSDs_%s_%s_AllSessions.json', patient_id_from_step1, neural_hemisphere);
output_params_path     = fullfile(preprocessed_output_folder, output_params_filename);

output_data_table = table();
if ~isempty(output_segments_list)
    try
        output_data_table = struct2table([output_segments_list{:}], 'AsArray', true);
    catch ME_struct2table
        error('struct2table failed: %s (segments=%d)', ME_struct2table.message, numel(output_segments_list));
    end
    fprintf('\nSaving CSV (%s)...\n', output_csv_filename);
    try
        writetable(output_data_table, output_csv_path);
        fprintf('CSV saved: %s (rows=%d)\n', output_csv_path, height(output_data_table));
    catch ME_writetable
        error('writetable failed: %s', ME_writetable.message);
    end
else
    warning('No segments generated. CSV not saved.');
end

% Save parameters JSON
fprintf('\nSaving parameters JSON...\n');
params_for_json.patient_id                           = patient_id_from_step1;
params_for_json.neural_hemisphere                    = neural_hemisphere;
params_for_json.pkg_data_hemisphere_used             = contralateral_pkg_hemisphere;
params_for_json.target_neural_segment_duration_sec   = target_neural_segment_duration_sec;
params_for_json.pkg_interpolation_interval_sec       = pkg_interpolation_interval_sec;
params_for_json.pwelch_config_used_in_matlab         = pwelch_config;
params_for_json.freq_ranges_defined_for_fooof        = freq_ranges_for_fooof;
params_for_json.preprocessing_steps_applied.high_pass_filter = apply_high_pass_filter;
if apply_high_pass_filter, params_for_json.preprocessing_details.high_pass_cutoff_Hz = high_pass_cutoff; end
params_for_json.preprocessing_steps_applied.stim_blanking   = use_stim_blanking;
params_for_json.source_step1_master_file             = step1_master_mat_filename;
params_for_json.source_pkg_file_used                 = selected_pkg_file;
params_for_json.utc_offset_hours_used                = UTCoffset_hours;
params_for_json.alignment_logic_description = sprintf('Center Alignment: Center of %.1fs neural segment aligned with Aligned_PKG_UnixTimestamp (Window=[T-%.1fs, T+%.1fs))', ...
    target_neural_segment_duration_sec, target_neural_segment_duration_sec/2, target_neural_segment_duration_sec/2);
params_for_json.script_name                          = mfilename('fullpath');
params_for_json.script_last_run_timestamp_utc        = datetime('now','TimeZone','UTC','Format','uuuu-MM-dd''T''HH:mm:ss.SSSZ');
params_for_json.matlab_version                       = version;
params_for_json.step2_output_folder                  = preprocessed_output_folder;
params_for_json.num_segments_generated               = height(output_data_table);

try
    json_string = jsonencode(params_for_json, 'PrettyPrint', true);
    fid = fopen(output_params_path, 'w');
    if fid == -1, error('Cannot create JSON file: %s', output_params_path); end
    fprintf(fid, '%s', json_string); fclose(fid);
    fprintf('JSON saved: %s\n', output_params_path);
catch ME_json
    warning('JSON save failed (%s). Saving MAT fallback.', ME_json.message);
    try
        save(strrep(output_params_path,'.json','.mat'), 'params_for_json');
    catch, warning('MAT fallback failed.'); end
end

fprintf('\nStep 2 completed (%.1fs neural windows, parallel over PKG epochs).\n', target_neural_segment_duration_sec);

%% ---- Helper Function: applyFilterToValidSegments_FIR ----
function [filtered_data, segments_too_short] = applyFilterToValidSegments_FIR(data, d, min_duration_sec_for_filter)
    % Applies FIR filter object d to non-NaN segments using filtfilt
    fs_filt = d.SampleRate;
    if ~isvector(data), error('Input data must be a vector.'); end
    if ~isnumeric(data) || ~isa(data,'double'), data = double(data); end

    nan_indices = isnan(data);
    transitions = diff([1; nan_indices(:); 1]);
    segment_starts = find(transitions == -1);
    segment_ends   = find(transitions == 1) - 1;

    filtered_data = data;
    segments_too_short = 0;

    min_samples_duration_based = round(min_duration_sec_for_filter * fs_filt);
    min_samples_for_filtfilt = 3 * d.FilterOrder;
    min_total_samples_required = max(min_samples_for_filtfilt, min_samples_duration_based);

    if isempty(segment_starts) && ~all(nan_indices)
        if length(data) < min_total_samples_required
            segments_too_short = 1; filtered_data(:) = NaN;
        else
            try, filtered_data = filtfilt(d, data); catch, filtered_data(:) = NaN; end
        end
        return;
    elseif all(nan_indices)
        return;
    end

    for i = 1:length(segment_starts)
        s = segment_starts(i); e = segment_ends(i);
        seg_len = e - s + 1;
        if seg_len < min_total_samples_required
            segments_too_short = segments_too_short + 1;
            filtered_data(s:e) = NaN; continue;
        end
        seg = data(s:e);
        try, filtered_data(s:e) = filtfilt(d, seg); catch, filtered_data(s:e) = NaN; end
    end
end
