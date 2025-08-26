%% Step2_github_parallel_clean.m  —  PARALLELIZED (parfor over PKG epochs)
%
% Preprocesses RC+S data from Step 1 Master, aligns CENTER of neural windows
% with PKG timestamps, computes PSDs, and streams results to CSV.
%
% Parallelization:
%   * parfor iterates over PKG epochs (thousands) instead of channels (few).
%   * Heavy read-only data shipped once per worker via parallel.pool.Constant.

clear; close all; clc;

%% 1) Configuration
fprintf('Starting Step 2\n');

script_run_folder = pwd;
project_base_path = '/home/jackson';

% --- Input Step 1 Master Data ---
target_patient_folder_step1  = 'RCS05L';
target_hemisphere_step1      = 'Left';
step1_master_mat_filename    = sprintf('Step1_MasterData_%s_%s_AllSessions_05.mat', ...
                                       target_patient_folder_step1, target_hemisphere_step1);
processed_data_folder_step1  = fullfile('/media/shortterm_ssd/jackson/RCS02L_Processed_Step1/');
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
    sprintf('step2_final', ...
    target_patient_folder_step1, target_hemisphere_step1));
if ~isfolder(preprocessed_output_folder), mkdir(preprocessed_output_folder); end
fprintf('Step 2 output will be saved to: %s\n', preprocessed_output_folder);

% --- PKG Data Parameters ---
pkg_interpolation_interval_sec = 30.0;
pkg_base_data_folder = fullfile('/media/shortterm_ssd/jackson/RCS05_(03JUL2019-05JUL2019)', 'PKG');
default_pkg_filename_pattern = 'scores_*.csv';

% --- Neural Preprocessing ---
apply_high_pass_filter = true;     % your original flag (set as you wish)
high_pass_cutoff       = 1;        % Hz
min_continuous_chunk_for_filter_sec = 5.0; %#ok<NASU>

% --- Segmentation and PSD Parameters ---
target_neural_segment_duration_sec = 30.0;   % CENTER-aligned 30 s window
pwelch_config.window_duration_sec  = 2.0;
pwelch_config.overlap_percent      = 50;
pwelch_config.nfft_multiplier      = 1;

% --- FOOOF ranges (metadata only here) ---
freq_ranges_for_fooof.LowFreq  = [10, 40];
freq_ranges_for_fooof.MidFreq  = [30, 90];
freq_ranges_for_fooof.WideFreq = [10, 90];

% --- Stim blanking flag (not used here) ---
use_stim_blanking = false;

% --- Initialize counters used later
num_rows_written = 0;

%% 2) Load Step1 Master Processed RC+S Data (no electrode_info required)
fprintf('\nLoading Step1 Master data from: %s\n', load_file_path_step1_master);
data_step1_master = load(load_file_path_step1_master);

required_master_vars = {'combinedDataTable_AllSessions','all_metaData_Master','target_patient_folder','target_hemisphere'};
missing_vars = setdiff(required_master_vars, fieldnames(data_step1_master));
if ~isempty(missing_vars)
    error('Missing required variable(s) "%s" from Step1 Master MAT file: %s.', strjoin(missing_vars, ', '), load_file_path_step1_master);
end

combinedDataTable      = data_step1_master.combinedDataTable_AllSessions;
all_metaData           = data_step1_master.all_metaData_Master;
patient_id_from_step1  = data_step1_master.target_patient_folder;
neural_hemisphere      = data_step1_master.target_hemisphere;

% Detect channels (support both legacy and new labels; optional 'contaxt' typo tolerated)
all_variable_names = combinedDataTable.Properties.VariableNames;
contact_channel_pattern = '^(?:key\d+_(?:contact|contaxt)_\d+_\d+|Contact_\d+_\d+)$';
channels_to_process = all_variable_names(~cellfun(@isempty, regexp(all_variable_names, contact_channel_pattern, 'once')));
fprintf('Detected %d contact channels to process.\n', numel(channels_to_process));

% Ensure Neural_UnixTimestamp exists (seconds)
if ~ismember('DerivedTime', combinedDataTable.Properties.VariableNames)
    error('Missing DerivedTime (ms) in master table.');
end
if ~ismember('Neural_UnixTimestamp', combinedDataTable.Properties.VariableNames)
    combinedDataTable.Neural_UnixTimestamp = combinedDataTable.DerivedTime / 1000;
    fprintf('  Created Neural_UnixTimestamp from DerivedTime.\n');
end
if ~isnumeric(combinedDataTable.Neural_UnixTimestamp)
    error('Neural_UnixTimestamp is not numeric.');
end

% UTC offset
if ~isempty(all_metaData) && iscell(all_metaData) && isfield(all_metaData{1}, 'UTCoffset')
    UTCoffset_hours = all_metaData{1}.UTCoffset;
elseif ~isempty(all_metaData) && isstruct(all_metaData) && numel(all_metaData)>=1 && isfield(all_metaData(1),'UTCoffset')
    UTCoffset_hours = all_metaData(1).UTCoffset;
else
    warning('UTCoffset not found; assuming 0.');
    UTCoffset_hours = 0;
end

%% 3) Load, Prepare, Interpolate PKG Data (Contralateral)
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

% Treat common empties as missing BEFORE typing to double
if ~isempty(valid_cols)
    opts = setvaropts(opts, valid_cols, 'TreatAsMissing', {'', 'NA', 'NaN', 'N/A', '.'});
    try, opts = setvaropts(opts, valid_cols, 'ThousandsSeparator', ','); catch, end
    for i = 1:numel(valid_cols)
        opts = setvartype(opts, valid_cols{i}, 'double');
    end
end

pkg_table_orig = readtable(pkg_file_path, opts);

% Post-read coercion guard
for i = 1:numel(valid_cols)
    c = valid_cols{i};
    if ~isnumeric(pkg_table_orig.(c))
        pkg_table_orig.(c) = str2double(string(pkg_table_orig.(c)));
    end
end

% Normalize BK sign as per your logic
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
            Aligned_PKG_UnixTimestamp = (t0:pkg_interpolation_interval_sec:t1)'; %#ok<NASGU>
            pkg_table_interp_aligned_scores = table(Aligned_PKG_UnixTimestamp);
            fprintf('Interpolating PKG scores:\n');
            for i = 1:numel(valid_cols)
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
                        pkg_table_interp_aligned_scores.(['Aligned_' nm]) = interp1(tt_u, yy_u, pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp, 'linear', nan);
                        fprintf('  - %s OK\n', nm);
                    else
                        pkg_table_interp_aligned_scores.(['Aligned_' nm]) = nan(height(pkg_table_interp_aligned_scores),1);
                        fprintf('  - %s skipped (not enough unique pts)\n', nm);
                    end
                else
                    pkg_table_interp_aligned_scores.(['Aligned_' nm]) = nan(height(pkg_table_interp_aligned_scores),1);
                    fprintf('  - %s skipped (not enough valid pts)\n', nm);
                end
            end
            if height(pkg_table_interp_aligned_scores)>0
                [~, idx_near] = min(abs(pkg_table_for_interp.PKG_UnixTimestamp_orig' - pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp), [], 2);
                pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str = pkg_table_for_interp.Date_Time(idx_near);
            end
        end
    end
end
fprintf('Interpolated PKG points: %d at %.1fs intervals.\n', height(pkg_table_interp_aligned_scores), pkg_interpolation_interval_sec);

%% 4) Neural & PKG time range sanity check
fprintf('\n🧠 Neural & PKG time range check:\n');
min_neural_time = min(combinedDataTable.Neural_UnixTimestamp);
max_neural_time = max(combinedDataTable.Neural_UnixTimestamp);
fprintf('  Neural: %.2f to %.2f (s)\n', min_neural_time, max_neural_time);
if ~isempty(pkg_table_interp_aligned_scores)
    min_pkg_time = min(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
    max_pkg_time = max(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
    fprintf('  PKG:    %.2f to %.2f (s)\n', min_pkg_time, max_pkg_time);
end

%% 4.5) Build read-only constants for parfor (PKG, time, channels, Welch)
% 1) Neural time (s) and fs estimate
t_neural = double(combinedDataTable.Neural_UnixTimestamp);
dt = diff(t_neural); dt = dt(isfinite(dt) & dt > 0);
if isempty(dt), error('Cannot estimate sampling rate from Neural_UnixTimestamp.'); end
fs_est = 1 / median(dt);

% 2) Welch parameters
win_samples = max(1, round(pwelch_config.window_duration_sec * fs_est));
noverlap   = round(pwelch_config.overlap_percent / 100 * win_samples);
nfft       = 2^nextpow2(max(win_samples, round(pwelch_config.nfft_multiplier * win_samples)));
w          = hann(win_samples, 'periodic');

welch_const       = parallel.pool.Constant(struct('w',w,'noverlap',noverlap,'nfft',nfft,'fs',fs_est,'win_samples',win_samples));
neural_time_const = parallel.pool.Constant(t_neural);

% 3) Channels (cell of vectors)
channels_cell = cell(numel(channels_to_process), 1);
for ch = 1:numel(channels_to_process)
    channels_cell{ch} = double(combinedDataTable.(channels_to_process{ch}));
end
chan_const = parallel.pool.Constant(channels_cell);

% 4) PKG struct
if ~isempty(pkg_table_interp_aligned_scores)
    pkg_struct = struct();
    pkg_struct.t  = double(pkg_table_interp_aligned_scores.Aligned_PKG_UnixTimestamp);
    if ismember('Aligned_PKG_DateTime_Str', pkg_table_interp_aligned_scores.Properties.VariableNames)
        pkg_struct.dt = pkg_table_interp_aligned_scores.Aligned_PKG_DateTime_Str;
    else
        pkg_struct.dt = repmat({''}, height(pkg_table_interp_aligned_scores), 1);
    end
    for i = 1:numel(valid_cols)
        nm = valid_cols{i};
        colnm = ['Aligned_' nm];
        if ismember(colnm, pkg_table_interp_aligned_scores.Properties.VariableNames)
            pkg_struct.(colnm) = pkg_table_interp_aligned_scores.(colnm);
        end
    end
else
    pkg_struct = struct('t', [], 'dt', {{}});
end
pkg_const = parallel.pool.Constant(pkg_struct);

fprintf('Parfor constants prepared: fs≈%.3f Hz, win=%d, noverlap=%d, nfft=%d, channels=%d, PKG epochs=%d.\n', ...
    fs_est, win_samples, noverlap, nfft, numel(channels_to_process), numel(pkg_struct.t));

%% 5) PARALLEL PROCESSING (capped workers + batched appends to CSV)
fprintf('\nProcessing neural PSDs with capped workers + batch appends...\n');

% Tunables
CAP_WORKERS   = 2;      % cap total workers
BATCH_EPOCHS  = 1500;    % PKG epochs per batch (trade RAM vs I/O)
tol_sec       = 1.5;     % tolerance around full window
max_nan_pct   = 90;      % skip segments too NaN-heavy

% Start (or cap) pool
try
    c = parcluster('local');
    nw = min(CAP_WORKERS, c.NumWorkers);
    if isempty(gcp('nocreate')), parpool(c, nw); end
catch ME
    warning('Could not start parpool (%s). Continuing in serial.', ME.message);
    nw = 1;
end
try, maxNumCompThreads(1); catch, end

% Prepare output CSV
output_csv_filename  = sprintf('Step2_Aligned120sPSDs_%s_%s_AllSessions.csv', patient_id_from_step1, neural_hemisphere);
output_csv_path      = fullfile(preprocessed_output_folder, output_csv_filename);
if isfile(output_csv_path), delete(output_csv_path); end

% Shortcuts for parfor
pkg       = pkg_const.Value;
nPKG      = numel(pkg.t);
half_win  = target_neural_segment_duration_sec / 2;
wc        = welch_const.Value;
t         = neural_time_const.Value;
chans     = chan_const.Value;
chns      = channels_to_process;

if nPKG == 0
    warning('No PKG-aligned timestamps inside neural range; nothing to process.');
end

first_write = true;

for s = 1:BATCH_EPOCHS:max(nPKG,1)
    e = min(s + BATCH_EPOCHS - 1, nPKG);
    thisN = max(e - s + 1, 0);

    seg_results = cell(thisN, 1);

    parfor k = 1:thisN
        ip = s + k - 1;
        center_t = pkg.t(ip);
        start_t  = center_t - half_win;
        end_t    = center_t + half_win;

        i1 = find(t >= start_t, 1, 'first');
        i2 = find(t <  end_t,  1, 'last');
        if isempty(i1) || isempty(i2) || i2 < i1
            seg_results{k} = {};
            continue;
        end

        seg_t = t(i1:i2);
        dur   = seg_t(end) - seg_t(1);
        if abs(dur - 2*half_win) > tol_sec
            seg_results{k} = {};
            continue;
        end

        entries = cell(1, numel(chns));
        for ch = 1:numel(chns)
            data = chans{ch};
            seg  = data(i1:i2);

            nan_pct = 100 * sum(isnan(seg)) / numel(seg);
            if nan_pct > max_nan_pct, continue; end

            vseg = seg(~isnan(seg));
            if numel(vseg) < wc.win_samples, continue; end

            try
                [Pxx, F] = pwelch(vseg, wc.w, wc.noverlap, wc.nfft, wc.fs);
            catch
                continue;
            end

            % Stringify (compact but precise)
            Pxx_str = strjoin(compose('%.8e', Pxx), ';');
            F_str   = strjoin(compose('%.4f', F),   ';');

            row = struct();
            row.PatientID  = patient_id_from_step1;
            row.Hemisphere = neural_hemisphere;
            row.Channel    = chns{ch};
            % Pretty label (normalize legacy to "Contact_A_B")
            row.ElectrodeLabel = regexprep(chns{ch}, '^key\d+_(?:contact|contaxt)_(\d+)_(\d+)$', 'Contact_$1_$2');

            row.Neural_Segment_Start_Unixtime = seg_t(1);
            row.Neural_Segment_End_Unixtime   = end_t;

            row.PSD_Data_Str         = Pxx_str;
            row.Frequency_Vector_Str = F_str;
            row.FS                   = wc.fs;

            row.Aligned_PKG_UnixTimestamp = center_t;
            if ~isempty(pkg.dt), row.Aligned_PKG_DateTime_Str = pkg.dt{ip}; else, row.Aligned_PKG_DateTime_Str = ''; end

            for sc = 1:numel(valid_cols)
                nm = ['Aligned_' valid_cols{sc}];
                if isfield(pkg, nm)
                    row.(nm) = pkg.(nm)(ip);
                end
            end

            entries{ch} = row;
        end

        entries = entries(~cellfun('isempty', entries));
        seg_results{k} = entries;
    end % parfor

    % ---- SERIAL: append this batch to CSV ----
    batch_list = [seg_results{:}];
    if ~isempty(batch_list)
        batch_table = struct2table([batch_list{:}], 'AsArray', true);

        try
            if first_write
                writetable(batch_table, output_csv_path);
                first_write = false;
            else
                writetable(batch_table, output_csv_path, 'WriteMode','append');
            end
        catch
            % Manual append for older MATLAB
            tmpf = [tempname '.csv'];
            writetable(batch_table, tmpf);
            fin  = fopen(tmpf,'r');
            if fin >= 0
                fout = fopen(output_csv_path,'a');
                if fout >= 0
                    if first_write
                        while true
                            tline = fgetl(fin);
                            if ~ischar(tline), break; end
                            fprintf(fout, '%s\n', tline);
                        end
                        first_write = false;
                    else
                        fgetl(fin); % skip header
                        while true
                            tline = fgetl(fin);
                            if ~ischar(tline), break; end
                            fprintf(fout, '%s\n', tline);
                        end
                    end
                    fclose(fout);
                else
                    warning('Could not open output CSV for append: %s', output_csv_path);
                end
                fclose(fin);
            else
                warning('Could not open temp file for append: %s', tmpf);
            end
            if exist(tmpf,'file'), delete(tmpf); end
        end

        num_rows_written = num_rows_written + height(batch_table);
    end

    clear seg_results batch_list batch_table
    if nPKG > 0
        fprintf('  Wrote batch %d–%d / %d (rows so far: %d)\n', s, e, nPKG, num_rows_written);
    end
end

fprintf('CSV saved: %s (rows=%d)\n', output_csv_path, num_rows_written);

% Expose count for JSON section
num_segments_generated = num_rows_written;

%% 6) Save parameters JSON
fprintf('\nSaving parameters JSON...\n');

output_params_filename = sprintf('Step2_params_120sPSDs_%s_%s_AllSessions.json', patient_id_from_step1, neural_hemisphere);
output_params_path     = fullfile(preprocessed_output_folder, output_params_filename);

params_for_json = struct();
params_for_json.patient_id                           = patient_id_from_step1;
params_for_json.neural_hemisphere                    = neural_hemisphere;
params_for_json.pkg_data_hemisphere_used             = contralateral_pkg_hemisphere;
params_for_json.target_neural_segment_duration_sec   = target_neural_segment_duration_sec;
params_for_json.pkg_interpolation_interval_sec       = pkg_interpolation_interval_sec;
params_for_json.pwelch_config_used_in_matlab         = pwelch_config;
params_for_json.freq_ranges_defined_for_fooof        = freq_ranges_for_fooof;
params_for_json.preprocessing_steps_applied          = struct();
params_for_json.preprocessing_steps_applied.high_pass_filter = apply_high_pass_filter;
if apply_high_pass_filter
    params_for_json.preprocessing_details = struct();
    params_for_json.preprocessing_details.high_pass_cutoff_Hz = high_pass_cutoff;
end
params_for_json.preprocessing_steps_applied.stim_blanking    = use_stim_blanking;
params_for_json.source_step1_master_file             = step1_master_mat_filename;
params_for_json.source_pkg_file_used                 = selected_pkg_file;
params_for_json.utc_offset_hours_used                = UTCoffset_hours;
params_for_json.alignment_logic_description          = sprintf('Center Alignment: Center of %.1fs neural segment aligned with Aligned_PKG_UnixTimestamp (Window=[T-%.1fs, T+%.1fs))', ...
    target_neural_segment_duration_sec, target_neural_segment_duration_sec/2, target_neural_segment_duration_sec/2);
params_for_json.script_name                          = mfilename('fullpath');
params_for_json.script_last_run_timestamp_utc        = datetime('now','TimeZone','UTC','Format','uuuu-MM-dd''T''HH:mm:ss.SSSZ');
params_for_json.matlab_version                       = version;
params_for_json.step2_output_folder                  = preprocessed_output_folder;
params_for_json.num_segments_generated               = num_segments_generated;

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
