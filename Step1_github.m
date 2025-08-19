%% Step 1: Process Raw RC+S Data (Revised for Multiple Sessions + Robust Samplerate)
% This script processes raw RC+S data from multiple session folders into MATLAB format,
% leveraging the Analysis-rcs-data toolbox. It generates individual MAT files
% for each session and a master MAT file with a combined data table.
% Includes robust handling for samplerate column data types before table concatenation.

clear all; close all; clc;

fprintf('Starting Step 1: Multi-Session RC+S Data Processing (with Robust Samplerate Handling)\n');
script_start_time = tic;

%% 1. Define paths and settings
% Set base path for project
base_path = fileparts(pwd);
%fprintf('Project base path set to: %s\n', base_path);

% --- Toolbox Path ---
toolbox_path = 'home/jackson/Analysis-rcs-data';
fprintf('Analysis-rcs-data toolbox found at: %s\n', toolbox_path);
addpath(genpath(toolbox_path)) 

%main_data_root = fullfile(base_path, 'RCS02 3 day sprint (05152019-05292019)/SCBS (Session1557950618572-1559185447897)');
main_data_root = '/media/shortterm_ssd/jackson/RCS05_(03JUL2019-05JUL2019)/SCBS/';
output_path = fullfile(base_path, 'jackson/step1_processed_data_multi_session_final_new/05');
if ~exist(output_path, 'dir')
    mkdir(output_path);
    fprintf('Created output directory: %s\n', output_path);
end
fprintf('Processed data will be saved to: %s\n', output_path);
vis_output_dir = fullfile(output_path, 'visualizations_per_session');
if ~exist(vis_output_dir, 'dir'), mkdir(vis_output_dir); end

%% 2. Define Target Patient and Hemisphere
target_patient_folder = 'RCS05L'; 
target_hemisphere = 'Left';   
process_rcs_flag = 3; 
default_TD_samplerate = 250; % Default Fs for TimeDomain if parsing fails
default_Accel_samplerate = 64; % Default Fs for Accelerometer if parsing fails (adjust as needed)
% Add other defaults if needed (e.g., Power, FFT)

patient_data_path = fullfile(main_data_root, target_hemisphere, target_patient_folder);
if ~exist(patient_data_path, 'dir'), error('Patient data path does not exist: %s', patient_data_path); end
fprintf('\nTargeting Patient: %s, Hemisphere: %s\nPath: %s\n', target_patient_folder, target_hemisphere, patient_data_path);

%% 3. Initialize Aggregators for Master File
combinedDataTable_AllSessions = table();
all_metaData_Master = {}; 
electrode_info_Master = struct(); 
electrode_info_extracted = false;
skipped_session_folders = {};
successfully_processed_sessions_details = table('Size',[0 3],'VariableTypes',{'string','string','string'},...
                                           'VariableNames',{'PatientFolder','SessionFolder','IndividualMatFilePath'});

%% 4. Identify and Loop Through Session Folders
session_folder_pattern = 'session*';
session_listing = dir(fullfile(patient_data_path, session_folder_pattern));
session_folders = session_listing([session_listing.isdir]);
%session_folders = session_folders(strcmp({session_folders.name}, ['session1558253403643']));

if isempty(session_folders)
    warning('No session folders found matching pattern "%s" in %s. Exiting.', session_folder_pattern, patient_data_path);
    return;
end
electrode_info_AllSessions = cell(1, length(session_folders));
contact_map = containers.Map();  % maps 'contact' → list of sessions/TD keys
fprintf('Found %d potential session folders for %s %s.\n', length(session_folders), target_patient_folder, target_hemisphere);

for i_sess = 1:length(session_folders)
    sessionName = session_folders(i_sess, :);
    current_session_name = session_folders(i_sess).name;
    current_session_path = fullfile(patient_data_path, current_session_name);
    individual_mat_filename = sprintf('Step1_SessData_%s_%s_%s.mat', target_patient_folder, target_hemisphere, current_session_name);
    individual_mat_filepath = fullfile(output_path, individual_mat_filename);


    fprintf('\n------------------------------------------------------\n');
    fprintf('Processing Session %d/%d: %s\n', i_sess, length(session_folders), current_session_name);

    clearvars unifiedDerivedTimes_sess timeDomainData_sess AccelData_sess PowerData_sess FFTData_sess AdaptiveData_sess ...
              timeDomainSettings_sess powerSettings_sess fftSettings_sess eventLogTable_sess metaData_sess ...
              stimSettingsOut_sess stimMetaData_sess stimLogSettings_sess DetectorSettings_sess ...
              AdaptiveStimSettings_sess AdaptiveEmbeddedRuns_StimSettings_sess combinedDataTable_sess;

    try
        device_npc_folder_listing = dir(fullfile(current_session_path, 'DeviceNPC*'));
        device_npc_folders = device_npc_folder_listing([device_npc_folder_listing.isdir]);
        if isempty(device_npc_folders)
            error('No DeviceNPC* folder found in %s.', current_session_path);
        elseif length(device_npc_folders) > 1
            warning('Multiple DeviceNPC* folders found in %s. Using the first one: %s.', current_session_path, device_npc_folders(1).name);
            device_raw_data_path = fullfile(current_session_path, device_npc_folders(1).name);
        else
            device_raw_data_path = fullfile(current_session_path, device_npc_folders(1).name);
        end
        fprintf('Device Raw Data Path: %s\n', device_raw_data_path);
        if ~exist(device_raw_data_path, 'dir'), error('Device folder does not exist: %s', device_raw_data_path); end
        json_files_check = dir(fullfile(device_raw_data_path, '*.json'));
        if isempty(json_files_check), error('No JSON files found in device folder: %s', device_raw_data_path); end
        

        fprintf('Running ProcessRCS for session %s...\n', current_session_name);
        [unifiedDerivedTimes_sess, timeDomainData_sess, ~, ~, ...
         AccelData_sess, ~, ~, PowerData_sess, ~, ~, FFTData_sess, ~, ~, AdaptiveData_sess, ~, ~, ...
         timeDomainSettings_sess, powerSettings_sess, fftSettings_sess, eventLogTable_sess, ...
         metaData_sess, stimSettingsOut_sess, stimMetaData_sess, stimLogSettings_sess, ...
         DetectorSettings_sess, AdaptiveStimSettings_sess, AdaptiveEmbeddedRuns_StimSettings_sess] = ProcessRCS(device_raw_data_path, process_rcs_flag);

        if isempty(timeDomainData_sess) && isempty(AccelData_sess) && isempty(PowerData_sess) 
            warning('ProcessRCS did not return significant data for session %s. Skipping.', current_session_name);
            skipped_session_folders{end+1} = current_session_name;
            continue; 
        end
        fprintf('ProcessRCS completed for session %s.\n', current_session_name);
        fprintf('  Available data streams for %s: TD:%s, Accel:%s, Power:%s, FFT:%s, Adaptive:%s\n', current_session_name, ...
            yesNoSymbol(~isempty(timeDomainData_sess)), yesNoSymbol(~isempty(AccelData_sess)), ...
            yesNoSymbol(~isempty(PowerData_sess)), yesNoSymbol(~isempty(FFTData_sess)), yesNoSymbol(~isempty(AdaptiveData_sess)));

        if ~isempty(unifiedDerivedTimes_sess)
            fprintf('Creating combinedDataTable for session %s...\n', current_session_name);
            electrode_info_sess = extract_electrode_channel_info(timeDomainSettings_sess, target_hemisphere, metaData_sess);
            electrode_info_AllSessions{i_sess} = electrode_info_sess;
            dataStreams_sess = {timeDomainData_sess, AccelData_sess, PowerData_sess, FFTData_sess, AdaptiveData_sess};
            combinedDataTable_sess = createCombinedTable(dataStreams_sess, unifiedDerivedTimes_sess, metaData_sess, electrode_info_sess);
            % === DEBUG: Check for all-NaN signal rows ===
            
            %contactCols = contains(combinedDataTable_sess.Properties.VariableNames, 'Contact_');
            contactCols = ~cellfun('isempty', regexp(combinedDataTable_sess.Properties.VariableNames, '^key\d+_contact_\d+_\d+$'));

            accelCols = contains(combinedDataTable_sess.Properties.VariableNames, 'Accel_');
            signalCols = contactCols | accelCols;
            
            nan_rows = all(ismissing(combinedDataTable_sess(:, signalCols)), 2);
            num_nan_rows = sum(nan_rows);
            
            % fprintf('⚠️  %d rows have all NaNs in TD and Accel signal columns for session %s.\n', ... 
            %      num_nan_rows, current_session_name);

            fprintf('Created combinedDataTable for session %s. Size: %d x %d\n', current_session_name, size(combinedDataTable_sess,1), size(combinedDataTable_sess,2));
            
            % --- Robust Samplerate Handling ---
            samplerate_cols_to_check = {
                'TD_samplerate', default_TD_samplerate; 
                'Accel_samplerate', default_Accel_samplerate
                % Add other samplerate columns here if they exist and cause issues, e.g.
                % 'Power_samplerate', default_Power_samplerate 
            };

            for sr_check_idx = 1:size(samplerate_cols_to_check, 1)
                col_name = samplerate_cols_to_check{sr_check_idx, 1};
                default_sr_val = samplerate_cols_to_check{sr_check_idx, 2};

                if ismember(col_name, combinedDataTable_sess.Properties.VariableNames)
                    if iscell(combinedDataTable_sess.(col_name))
                        fprintf('  Session %s: Column "%s" is a cell. Attempting conversion to numeric.\n', current_session_name, col_name);
                        numeric_sr_col = NaN(height(combinedDataTable_sess), 1);
                        for sr_row_idx = 1:height(combinedDataTable_sess)
                            cell_content = combinedDataTable_sess.(col_name){sr_row_idx};
                            if ischar(cell_content) || isstring(cell_content)
                                val = str2double(regexp(cell_content, '\d+\.?\d*', 'match', 'once'));
                                if ~isempty(val) && ~isnan(val)
                                    numeric_sr_col(sr_row_idx) = val;
                                else
                                    numeric_sr_col(sr_row_idx) = default_sr_val;
                                    fprintf('    Warning: Could not parse SR from "%s" for %s in session %s, row %d. Defaulting to %.0f Hz.\n', char(cell_content), col_name, current_session_name, sr_row_idx, default_sr_val);
                                end
                            elseif isnumeric(cell_content)
                                numeric_sr_col(sr_row_idx) = cell_content;
                            else
                                numeric_sr_col(sr_row_idx) = default_sr_val;
                                fprintf('    Warning: Unexpected cell content for %s in session %s, row %d. Defaulting to %.0f Hz.\n', col_name, current_session_name, sr_row_idx, default_sr_val);
                            end
                        end
                        combinedDataTable_sess.(col_name) = numeric_sr_col;
                        fprintf('  Session %s: Column "%s" converted to numeric.\n', current_session_name, col_name);
                    elseif ~isnumeric(combinedDataTable_sess.(col_name))
                        fprintf('  Session %s: Column "%s" is not numeric and not cell. Attempting direct conversion or defaulting.\n', current_session_name, col_name);
                        try
                            combinedDataTable_sess.(col_name) = double(combinedDataTable_sess.(col_name));
                        catch
                            fprintf('    Warning: Could not convert %s to double for session %s. Defaulting column to %.0f Hz.\n', col_name, current_session_name, default_sr_val);
                            combinedDataTable_sess.(col_name) = repmat(default_sr_val, height(combinedDataTable_sess), 1);
                        end
                    end
                    if any(isnan(combinedDataTable_sess.(col_name)))
                        nan_indices_sr = isnan(combinedDataTable_sess.(col_name));
                        combinedDataTable_sess.(col_name)(nan_indices_sr) = default_sr_val;
                        fprintf('  Session %s: NaN values in %s defaulted to %.0f Hz.\n', current_session_name, col_name, default_sr_val);
                    end
                else % Column does not exist, add it with default if its corresponding data stream likely exists
                    % Example: If Accel_XSamples exists, Accel_samplerate should ideally too.
                    if strcmp(col_name, 'Accel_samplerate') && ismember('Accel_XSamples', combinedDataTable_sess.Properties.VariableNames)
                         fprintf('  Warning: Column %s not found in combinedDataTable_sess for session %s. Adding with default %.0f Hz.\n', col_name, current_session_name, default_sr_val);
                         combinedDataTable_sess.(col_name) = repmat(default_sr_val, height(combinedDataTable_sess), 1);
                    elseif strcmp(col_name, 'TD_samplerate') && ismember('TD_key0', combinedDataTable_sess.Properties.VariableNames) % If TD data exists
                         fprintf('  Warning: Column %s not found for session %s. Adding with default %.0f Hz.\n', col_name, current_session_name, default_sr_val);
                         combinedDataTable_sess.(col_name) = repmat(default_sr_val, height(combinedDataTable_sess), 1);
                    end
                end
            end

        else
            warning('unifiedDerivedTimes is empty for session %s. Cannot create combinedDataTable_sess. Skipping session.', current_session_name);
            skipped_session_folders{end+1} = current_session_name;
            continue;
        end

        % if ~electrode_info_extracted && ~isempty(timeDomainSettings_sess)
        %     fprintf('Extracting electrode channel information from session %s (first success)...\n', current_session_name);
        %     electrode_info_Master = extract_electrode_channel_info(timeDomainSettings_sess, target_hemisphere, metaData_sess); 
        %     if ~isempty(fieldnames(electrode_info_Master))
        %         electrode_info_extracted = true;
        %         fprintf('Electrode information extracted and stored.\n');
        %     else
        %         fprintf('Electrode information extraction yielded no channels. Will try next session.\n');
        %     end
        % end
        if ~isempty(timeDomainSettings_sess)
            fprintf('Extracting electrode channel information from session %s...\n', current_session_name);
            
            fprintf('  Contact configuration for session %s:\n', current_session_name);
            td_keys = fieldnames(electrode_info_sess);
            for k = 1:length(td_keys)
                chan = td_keys{k};
                contact = electrode_info_sess.(chan).contacts;
                fs = electrode_info_sess.(chan).fs;
                fprintf('    %s → Contact: %s | Fs: %s\n', chan, contact, num2str(fs));
            end
            % 🔁 Add contact entries to map
            td_keys = fieldnames(electrode_info_sess);
            for k = 1:length(td_keys)
                tdkey = td_keys{k};
                contact = electrode_info_sess.(tdkey).contacts;
                fs = electrode_info_sess.(tdkey).fs;
                if strcmpi(contact, 'Disabled') || strcmpi(contact, 'Unknown')
                    continue;
                end

                if ~isKey(contact_map, contact)
                    contact_map(contact) = {};
                end
                contact_entries = contact_map(contact);
                
     
                contact_entries{end+1} = struct( ...
                    'session_idx', i_sess, ...
                    'td_key', tdkey, ...
                    'fs', fs, ...
                    'key_idx', electrode_info_sess.(tdkey).key_idx ...
                );
                contact_map(contact) = contact_entries;  % ← reassign the updated cell array
            end
        end
        % if ~isempty(timeDomainData_sess)
        %     fprintf('Creating visualization for session %s...\n', current_session_name);
        %     try
        %         create_broadband_visualization(timeDomainData_sess, electrode_info_sess, vis_output_dir, current_session_name, device_raw_data_path, target_patient_folder); 
        %     catch vis_ME
        %          warning('Failed to create visualization for session %s: %s', current_session_name, vis_ME.message);
        %     end
        % end

        if ~isempty(combinedDataTable_sess)
            if isempty(combinedDataTable_AllSessions.Properties.VariableNames)
                combinedDataTable_AllSessions = combinedDataTable_sess;
            else
                current_cols = combinedDataTable_sess.Properties.VariableNames;
                master_cols = combinedDataTable_AllSessions.Properties.VariableNames;
                
                % Align columns: Add missing columns from master to current, and from current to master, fill with appropriate NaNs/empty cells
                % This is a more robust way to prepare for vertcat if column sets differ slightly.
                
                % Add missing columns from master to current session's table
                missing_in_current = setdiff(master_cols, current_cols);
                for mc = 1:length(missing_in_current)
                    col_to_add = missing_in_current{mc};
                    % Determine type from master table for appropriate fill
                    if isnumeric(combinedDataTable_AllSessions.(col_to_add)) || islogical(combinedDataTable_AllSessions.(col_to_add))
                        combinedDataTable_sess.(col_to_add) = NaN(height(combinedDataTable_sess), size(combinedDataTable_AllSessions.(col_to_add),2) );
                    elseif iscell(combinedDataTable_AllSessions.(col_to_add))
                        combinedDataTable_sess.(col_to_add) = cell(height(combinedDataTable_sess), size(combinedDataTable_AllSessions.(col_to_add),2) );
                    elseif isdatetime(combinedDataTable_AllSessions.(col_to_add))
                        combinedDataTable_sess.(col_to_add) = NaT(height(combinedDataTable_sess), size(combinedDataTable_AllSessions.(col_to_add),2) );
                    else % Default to NaN for other types or handle specifically
                        combinedDataTable_sess.(col_to_add) = NaN(height(combinedDataTable_sess), size(combinedDataTable_AllSessions.(col_to_add),2) );
                    end
                end

                % Add missing columns from current to master session's table
                missing_in_master = setdiff(current_cols, master_cols);
                for mc = 1:length(missing_in_master)
                    col_to_add = missing_in_master{mc};
                     if isnumeric(combinedDataTable_sess.(col_to_add)) || islogical(combinedDataTable_sess.(col_to_add))
                        combinedDataTable_AllSessions.(col_to_add) = NaN(height(combinedDataTable_AllSessions), size(combinedDataTable_sess.(col_to_add),2) );
                    elseif iscell(combinedDataTable_sess.(col_to_add))
                        combinedDataTable_AllSessions.(col_to_add) = cell(height(combinedDataTable_AllSessions), size(combinedDataTable_sess.(col_to_add),2) );
                    elseif isdatetime(combinedDataTable_sess.(col_to_add))
                        combinedDataTable_AllSessions.(col_to_add) = NaT(height(combinedDataTable_AllSessions), size(combinedDataTable_sess.(col_to_add),2) );
                    else
                        combinedDataTable_AllSessions.(col_to_add) = NaN(height(combinedDataTable_AllSessions), size(combinedDataTable_sess.(col_to_add),2) );
                    end
                end
                
                % Now ensure order is the same before vertcat
                if ~isempty(combinedDataTable_AllSessions.Properties.VariableNames) % if master wasn't initially empty
                    combinedDataTable_sess = combinedDataTable_sess(:, combinedDataTable_AllSessions.Properties.VariableNames);
                elseif ~isempty(combinedDataTable_sess.Properties.VariableNames) % master was empty, current is not, set master to current structure
                     combinedDataTable_AllSessions = combinedDataTable_sess([],:); % Empty table with same columns as current
                     combinedDataTable_AllSessions = combinedDataTable_AllSessions(:, combinedDataTable_sess.Properties.VariableNames); % ensure order
                end


                try
                    combinedDataTable_AllSessions = vertcat(combinedDataTable_AllSessions, combinedDataTable_sess);
                    fprintf('Appended combinedDataTable from session %s to master table.\n', current_session_name);
                catch vertcat_ME
                    warning('VERTCAT FAILED for session %s despite type and column checks. Error: %s. Skipping append.', current_session_name, vertcat_ME.message);
                    fprintf('Master table cols: %s\n', strjoin(combinedDataTable_AllSessions.Properties.VariableNames,', '));
                    fprintf('Session table cols: %s\n', strjoin(combinedDataTable_sess.Properties.VariableNames,', '));
                    skipped_session_folders{end+1} = current_session_name; % Add to skipped
                    continue; % Skip to next session
                end
            end
        end
        
        metaData_sess.sessionName = current_session_name; 
        all_metaData_Master{end+1} = metaData_sess;

        individual_mat_filename = sprintf('Step1_SessData_%s_%s_%s.mat', target_patient_folder, target_hemisphere, current_session_name);
        individual_mat_filepath = fullfile(output_path, individual_mat_filename);


        fprintf('Saving individual MAT file for session %s to: %s\n', current_session_name, individual_mat_filepath);
        save(individual_mat_filepath, ...
            'unifiedDerivedTimes_sess', 'timeDomainData_sess', 'AccelData_sess', 'PowerData_sess', 'FFTData_sess', 'AdaptiveData_sess', ...
            'timeDomainSettings_sess', 'powerSettings_sess', 'fftSettings_sess', 'eventLogTable_sess', 'metaData_sess', ...
            'stimSettingsOut_sess', 'stimMetaData_sess', 'stimLogSettings_sess', 'DetectorSettings_sess', ...
            'AdaptiveStimSettings_sess', 'AdaptiveEmbeddedRuns_StimSettings_sess', ...
            'combinedDataTable_sess', ... 
            'current_session_name', 'target_patient_folder', 'target_hemisphere', ...
             '-v7.3');
        fprintf('Successfully saved individual MAT for session %s.\n', current_session_name);
        
        new_row = {string(target_patient_folder), string(current_session_name), string(individual_mat_filepath)};
        successfully_processed_sessions_details = [successfully_processed_sessions_details; new_row];

    catch ME
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('ERROR processing session %s: %s\n', current_session_name, ME.message);
        fprintf('Error details: %s\n', ME.getReport('extended', 'hyperlinks','off'));
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        skipped_session_folders{end+1} = current_session_name;
    end 
end 



%% Step 1.5: Rebuild Master MAT from Individual Session Files (Parallel, no electrode_info)

clear; clc;
fprintf('Rebuilding master MAT file from individual sessions (parallel)...\n');

% --- Configuration ---
output_path            = '/home/jackson/step1_processed_data_multi_session_final_new/05'; % Update as needed
target_patient_folder  = 'RCS05R';
target_hemisphere      = 'Right';

% Find all individual .mat files
mat_files = dir(fullfile(output_path, sprintf('Step1_SessData_%s_%s_*.mat', target_patient_folder, target_hemisphere)));
if isempty(mat_files), error('No individual session .mat files found in: %s', output_path); end
nFiles = numel(mat_files);

% --- Start/ensure a local pool
try
    pool = gcp('nocreate');
    if isempty(pool), parpool('local'); end
catch ME
    warning('Could not start a parallel pool: %s\nProceeding single-threaded.', ME.message);
end

% --- Containers for parallel import
sess_ok          = false(nFiles,1);
sess_err_msg     = strings(nFiles,1);
sess_paths       = strings(nFiles,1);
sess_names       = strings(nFiles,1);

T_sess_all       = cell(nFiles,1);   % combinedDataTable_sess
meta_sess_all    = cell(nFiles,1);   % metaData_sess (optional)
varnames_all     = cell(nFiles,1);   % varnames for union

% --- Parallel load of all sessions
parfor i = 1:nFiles
    name_i = mat_files(i).name;
    path_i = fullfile(mat_files(i).folder, name_i);

    local_ok  = false;
    local_err = "";

    try
        S = load(path_i);

        if ~isfield(S,'combinedDataTable_sess') || isempty(S.combinedDataTable_sess)
            error('No combinedDataTable_sess');
        end

        % Session name from filename
        match = regexp(name_i, 'Step1_SessData_.*?_(session\d+)\.mat', 'tokens');
        sessionFolderName = '';
        if ~isempty(match), sessionFolderName = match{1}{1}; end

        % Populate outputs
        T_sess_all{i}   = S.combinedDataTable_sess;
        varnames_all{i} = S.combinedDataTable_sess.Properties.VariableNames;

        if isfield(S,'metaData_sess'), meta_sess_all{i} = S.metaData_sess; end

        sess_paths(i) = string(path_i);
        sess_names(i) = string(sessionFolderName);

        local_ok = true;
    catch ME
        local_err = string(ME.message);
    end

    sess_ok(i)      = local_ok;
    sess_err_msg(i) = local_err;
end

% Report failures (but continue with successes)
if any(~sess_ok)
    bad_idx = find(~sess_ok);
    fprintf('⚠️ %d/%d sessions failed to load:\n', numel(bad_idx), nFiles);
    for b = bad_idx(:).'
        fprintf('   - %s : %s\n', mat_files(b).name, sess_err_msg(b));
    end
end

good_idx = find(sess_ok);
if isempty(good_idx)
    error('All sessions failed to load. Aborting.');
end
fprintf('✅ %d/%d sessions loaded successfully.\n', numel(good_idx), nFiles);

% --- Initialize master containers
combinedDataTable_AllSessions = table();
all_metaData_Master           = {};
successfully_processed_sessions_details = table('Size', [0 3], ...
    'VariableTypes', {'string','string','string'}, ...
    'VariableNames', {'PatientFolder','SessionFolder','IndividualMatFilePath'});

% --- Build union of variable names across successful sessions
union_vars = {};
seen = containers.Map('KeyType','char','ValueType','logical');
for ii = good_idx(:).'
    vns = varnames_all{ii};
    for v = 1:numel(vns)
        vn = vns{v};
        if ~isKey(seen, vn)
            union_vars{end+1} = vn; %#ok<AGROW>
            seen(vn) = true;
        end
    end
end

% --- Merge sessions serially (align to union of columns)
for ii = good_idx(:).'
    T = T_sess_all{ii};

    % Add missing columns to this session table (types inferred from another session that has it)
    missing_cols = setdiff(union_vars, T.Properties.VariableNames);
    for c = 1:numel(missing_cols)
        col = missing_cols{c};
        % Find an example column from another session to infer type
        example_val = [];
        for jj = good_idx(:).'
            if ismember(col, varnames_all{jj})
                example_val = T_sess_all{jj}.(col);
                break;
            end
        end
        nRows = height(T);
        if isempty(example_val)
            % Fallback if no example found
            T.(col) = strings(nRows,1);
        else
            if isnumeric(example_val)
                T.(col) = NaN(nRows,1);
            elseif iscell(example_val)
                T.(col) = cell(nRows,1);
            elseif isdatetime(example_val)
                T.(col) = NaT(nRows,1);
            elseif islogical(example_val)
                T.(col) = false(nRows,1);
            elseif isstring(example_val) || ischar(example_val)
                T.(col) = strings(nRows,1);
            else
                T.(col) = strings(nRows,1);
            end
        end
    end

    % Ensure consistent column order
    T = T(:, union_vars);

    % Append
    if isempty(combinedDataTable_AllSessions)
        combinedDataTable_AllSessions = T;
    else
        combinedDataTable_AllSessions = [combinedDataTable_AllSessions; T]; %#ok<AGROW>
    end

    % Collect metadata
    if ~isempty(meta_sess_all{ii})
        all_metaData_Master{end+1} = meta_sess_all{ii}; %#ok<AGROW>
    end

    % Log success table row
    new_row = {string(target_patient_folder), sess_names(ii), sess_paths(ii)};
    successfully_processed_sessions_details = [successfully_processed_sessions_details; new_row]; %#ok<AGROW>
end

% --- Optional: sort by time if available
if ismember('localTime', combinedDataTable_AllSessions.Properties.VariableNames)
    try
        combinedDataTable_AllSessions = sortrows(combinedDataTable_AllSessions, 'localTime');
    catch
        % If localTime is heterogeneous type across sessions, skip sorting
        warning('Could not sort by localTime due to mixed types.');
    end
end

% --- Save rebuilt master .mat (NO electrode_info variables)
master_mat_filename = sprintf('Step1_MasterData_%s_%s_AllSessions_05.mat', target_patient_folder, target_hemisphere);
master_mat_filepath = fullfile(output_path, master_mat_filename);
save(master_mat_filepath, ...
    'combinedDataTable_AllSessions', 'all_metaData_Master', ...
    'target_patient_folder', 'target_hemisphere', ...
    'successfully_processed_sessions_details', ...
    'output_path', ...
    '-v7.3');
fprintf('✅ Rebuilt master MAT file saved to:\n%s\n', master_mat_filepath);

%% Helper Function: yesNoSymbol 
function symbol = yesNoSymbol(condition)
    if condition, symbol = '✓'; else, symbol = '✗'; end
end

function [combinedDataTable] = createCombinedTable(dataStreams, unifiedDerivedTimes, metaData, electrode_info)
    % --- Always use full unifiedDerivedTimes as the time base ---
    combinedDataTable = table(unifiedDerivedTimes, 'VariableNames', {'DerivedTime'});
    numRows = height(combinedDataTable);

    % TD signal mapping (padded)
    if ~isempty(dataStreams{1}) && ismember('key0', dataStreams{1}.Properties.VariableNames)
        tdData = dataStreams{1};
        for iKey = 0:3
            keyName = sprintf('key%d', iKey); 
            chan_name = sprintf('TD_key%d', iKey);
            if ~isfield(electrode_info, chan_name), continue; end
            contact_label = electrode_info.(chan_name).contacts;
            if strcmpi(contact_label, 'Disabled') || strcmpi(contact_label, 'Unknown'), continue; end
            %contact_varname = matlab.lang.makeValidName(sprintf('Contact_%s', strrep(contact_label, '-', '_')));
            contact_varname = matlab.lang.makeValidName(sprintf('%s_contact_%s', keyName, strrep(contact_label, '-', '_')));


            if ismember(keyName, tdData.Properties.VariableNames)
                [~, ia, ib] = intersect(unifiedDerivedTimes, tdData.DerivedTime, 'stable');
                padded_signal = NaN(numRows, 1);
                padded_signal(ia) = tdData.(keyName)(ib);
                combinedDataTable.(contact_varname) = padded_signal;
                fprintf('  ✅ Assigned and padded %s → %s (mapped %d/%d samples)\n', ...
                    keyName, contact_varname, length(ib), numRows);
                % Pad samplerate if present
                if ismember('samplerate', tdData.Properties.VariableNames)
                    padded_sr = NaN(numRows, 1);
                    padded_sr(ia) = tdData.samplerate(ib);
                    combinedDataTable.TD_samplerate = padded_sr;
                else
                    combinedDataTable.TD_samplerate = repmat(500, numRows, 1); % fallback default
                end
            end
        end
    end

    % --- Accel Data Padding ---
    accelData = dataStreams{2};
    if ~isempty(accelData) && ismember('Accel_XSamples', accelData.Properties.VariableNames) && ismember('DerivedTime', accelData.Properties.VariableNames)
        [~, ia, ib] = intersect(unifiedDerivedTimes, accelData.DerivedTime, 'stable');
        % X/Y/Z axes
        accelFields = {'Accel_XSamples', 'Accel_YSamples', 'Accel_ZSamples'};
        for f = 1:length(accelFields)
            fld = accelFields{f};
            padded_accel = NaN(numRows, 1);
            if ismember(fld, accelData.Properties.VariableNames)
                padded_accel(ia) = accelData.(fld)(ib);
            end
            combinedDataTable.(fld) = padded_accel;
        end
    
        % Accel_samplerate logic
        if ismember('samplerate', accelData.Properties.VariableNames)
            sr_vec = accelData.samplerate;
            if iscell(sr_vec), sr_vec = cellfun(@double, sr_vec); end
            sr_unique = unique(sr_vec(~isnan(sr_vec)));
            if isempty(sr_unique)
                session_sr = 64;
                fprintf('⚠️  No valid Accel_samplerate found, using default 64 Hz\n');
            else
                session_sr = sr_unique(1);
                if numel(sr_unique) > 1
                    fprintf('⚠️  Multiple Accel_samplerate values detected ([%s]) — using first: %.0f Hz\n', num2str(sr_unique'), session_sr);
                end
            end
            combinedDataTable.Accel_samplerate = repmat(session_sr, numRows, 1);
        else
            combinedDataTable.Accel_samplerate = repmat(64, numRows, 1); % fallback default
            fprintf('⚠️  No Accel_samplerate column, defaulting to 64 Hz\n');
        end
    end

    % Add localTime conversion
    if ~isempty(combinedDataTable) && ismember('DerivedTime', combinedDataTable.Properties.VariableNames) && ~isempty(metaData) && isfield(metaData, 'UTCoffset')
        timeFormat = sprintf('%+03.0f:00',metaData.UTCoffset);
        if isnumeric(combinedDataTable.DerivedTime)
            localTime = datetime(combinedDataTable.DerivedTime/1000, 'ConvertFrom','posixTime','TimeZone',timeFormat,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            combinedDataTable = addvars(combinedDataTable,localTime,'Before',1);
        else
            warning('DerivedTime non-numeric in createCombinedTable.'); 
        end
    else
        if isempty(combinedDataTable) || ~ismember('DerivedTime', combinedDataTable.Properties.VariableNames)
            warning('Cannot add localTime: DerivedTime missing/empty.');
        else
            warning('Cannot add localTime: metaData/UTCoffset missing.'); 
        end
    end

    fprintf('\n✅ Combined table created: %d rows x %d columns\n', height(combinedDataTable), width(combinedDataTable));
end


%% Helper Function: extract_electrode_channel_info 
function electrode_info = extract_electrode_channel_info(timeDomainSettings, hemisphere, metaData)
    fprintf('Extracting electrode channel information...\n');
    electrode_info = struct();
    if isempty(timeDomainSettings) || ~istable(timeDomainSettings) || height(timeDomainSettings) == 0
        fprintf('  No time domain settings available. Cannot extract electrode info.\n'); return;
    end
    if ismember('TDsettings', timeDomainSettings.Properties.VariableNames)
        firstValidSettingsRow = find(~cellfun(@isempty, timeDomainSettings.TDsettings), 1, 'first');
        if ~isempty(firstValidSettingsRow)
            tdSettingsArray = timeDomainSettings.TDsettings{firstValidSettingsRow};
            fprintf('  Found TDsettings with %d configurations.\n', length(tdSettingsArray));
            for chan_idx = 1:length(tdSettingsArray) 
                chan_config = tdSettingsArray(chan_idx); 
                chan_name = sprintf('TD_key%d', chan_idx-1);
                electrode_info.(chan_name) = struct('type','TD','hemisphere',hemisphere,'fs',NaN,'contacts','Unknown','location','Unknown');
                electrode_info.(chan_name).key_idx = chan_idx - 1;
                if isfield(chan_config, 'sampleRate') && ~isempty(chan_config.sampleRate) && ...
                   ~strcmpi(chan_config.sampleRate, 'disabled') && ~strcmpi(chan_config.sampleRate, 'unexpected')
                    if isnumeric(chan_config.sampleRate), fs_val = chan_config.sampleRate;
                    else, fs_val_match = regexp(chan_config.sampleRate, '\d+', 'match', 'once');
                        if ~isempty(fs_val_match), fs_val = str2double(fs_val_match); else, fs_val = NaN; end
                    end
                    electrode_info.(chan_name).fs = fs_val;
                    if isfield(chan_config, 'plusInput') && isfield(chan_config, 'minusInput')
                        electrode_info.(chan_name).contacts = sprintf('%s-%s', chan_config.plusInput, chan_config.minusInput);
                        plusNumStr = regexp(chan_config.plusInput, '\d+', 'match', 'once');
                        minusNumStr = regexp(chan_config.minusInput, '\d+', 'match', 'once');
                        if ~isempty(plusNumStr) && ~isempty(minusNumStr)
                            plusNum = str2double(plusNumStr); minusNum = str2double(minusNumStr);
                            if max(plusNum, minusNum) <= 3, electrode_info.(chan_name).location = 'STN/DBS';
                            elseif min(plusNum, minusNum) >= 4, electrode_info.(chan_name).location = 'Cortical/ECoG';
                            else, electrode_info.(chan_name).location = 'Mixed/Unknown'; end
                        end
                    end
                else, electrode_info.(chan_name).contacts = 'Disabled'; electrode_info.(chan_name).location = 'Disabled'; end
            end
        else, fprintf('  TDsettings field present but all entries are empty.\n'); end
    elseif any(ismember({'chan1','chan2','chan3','chan4'}, timeDomainSettings.Properties.VariableNames))
        fprintf('  Using alternative approach (chanX fields) for electrode info.\n');
        for chan_num = 1:4 
            chan_field_name = sprintf('chan%d', chan_num); td_key_name = sprintf('TD_key%d', chan_num-1);
            electrode_info.(td_key_name) = struct('type','TD','hemisphere',hemisphere,'fs',NaN,'contacts','Unknown','location','Unknown');
            if ismember(chan_field_name, timeDomainSettings.Properties.VariableNames)
                config_str = ''; fs_val_row = NaN;
                for row = 1:height(timeDomainSettings)
                    if iscell(timeDomainSettings.(chan_field_name)) && ~isempty(timeDomainSettings.(chan_field_name){row})
                        config_str = timeDomainSettings.(chan_field_name){row};
                        if ismember('samplingRate', timeDomainSettings.Properties.VariableNames) && ~isnan(timeDomainSettings.samplingRate(row))
                             fs_val_row = timeDomainSettings.samplingRate(row); break; 
                        end
                    end
                end
                if ~isnan(fs_val_row), electrode_info.(td_key_name).fs = fs_val_row; end
                if ~isempty(config_str)
                    contacts_match = regexp(config_str, '\+(\w+)-(\w+)', 'tokens', 'once');
                    if ~isempty(contacts_match)
                        electrode_info.(td_key_name).contacts = sprintf('%s-%s', contacts_match{1}, contacts_match{2});
                        % (Location logic as above)
                    end
                end
            end
        end
    else, fprintf('  Neither "TDsettings" nor "chanX" fields found for electrode info.\n'); end
    % (Display logic as before)
end

%% Helper Function: create_broadband_visualization 
function create_broadband_visualization(timeDomainData_sess, electrode_info_Master, vis_dir, session_name, device_path_sess, patient_folder)
    % (Code from previous response, ensure it's correct and handles missing electrode_info_Master gracefully)
    fprintf('  Attempting to create broadband signal visualization...\n');
    if isempty(timeDomainData_sess) || height(timeDomainData_sess) < 2, fprintf('  No TD data for vis.\n'); return; end
    if ~ismember('DerivedTime', timeDomainData_sess.Properties.VariableNames), fprintf('  DerivedTime missing for vis.\n'); return; end
    
    total_duration_sec = (timeDomainData_sess.DerivedTime(end) - timeDomainData_sess.DerivedTime(1)) / 1000;
    vis_duration_sec = min(600, total_duration_sec); 
    
    channel_cols_to_plot = {}; td_keys = {};
    if ~isempty(fieldnames(electrode_info_Master)),  td_keys = fieldnames(electrode_info_Master); end

    for k_idx = 1:length(td_keys)
        key_name = td_keys{k_idx}; 
        if isfield(electrode_info_Master.(key_name), 'fs') && ~isnan(electrode_info_Master.(key_name).fs) 
             actual_col_name = strrep(key_name, 'TD_', ''); 
             if ismember(actual_col_name, timeDomainData_sess.Properties.VariableNames)
                channel_cols_to_plot{end+1} = actual_col_name;
             end
        end
    end
    if isempty(channel_cols_to_plot) % Fallback if electrode_info_Master is empty or no channels matched
        for i=0:3, potential_col = sprintf('key%d',i); if ismember(potential_col, timeDomainData_sess.Properties.VariableNames), channel_cols_to_plot{end+1}=potential_col; end; end
    end
    if isempty(channel_cols_to_plot), fprintf('  No channels found for plotting.\n'); return; end
    
    time_end_idx = find(timeDomainData_sess.DerivedTime >= (timeDomainData_sess.DerivedTime(1) + vis_duration_sec*1000), 1, 'first');
    if isempty(time_end_idx), time_end_idx = height(timeDomainData_sess); end
    
    max_points_vis = 20000; num_points_in_window = time_end_idx; downsample_factor_vis = 1;
    if num_points_in_window > max_points_vis, downsample_factor_vis = ceil(num_points_in_window / max_points_vis); end
    idx_to_plot = 1:downsample_factor_vis:time_end_idx;
    
    fig_vis = figure('Position', [100, 100, 1200, max(400, 200*length(channel_cols_to_plot))], 'Visible', 'off'); 
    for i_ch_plot = 1:length(channel_cols_to_plot)
        subplot(length(channel_cols_to_plot), 1, i_ch_plot);
        current_col_name = channel_cols_to_plot{i_ch_plot}; 
        y_data = timeDomainData_sess.(current_col_name)(idx_to_plot);
        time_data_ms = timeDomainData_sess.DerivedTime(idx_to_plot);
        valid_idx_plot = ~isnan(y_data);
        if sum(valid_idx_plot) > 1, plot(time_data_ms(valid_idx_plot)/1000, y_data(valid_idx_plot)); 
        else, text(0.5, 0.5, 'No valid data', 'HorizontalAlignment', 'center'); end
        
        chan_key_for_title = ['TD_' current_col_name]; 
        title_str = sprintf('Channel %s', current_col_name);
        if ~isempty(fieldnames(electrode_info_Master)) && isfield(electrode_info_Master, chan_key_for_title)
            info = electrode_info_Master.(chan_key_for_title);
            if isfield(info, 'location') && ~strcmp(info.location, 'Unknown'), title_str = [title_str ' - ' info.location]; end
            if isfield(info, 'contacts') && ~strcmp(info.contacts, 'Unknown'), title_str = [title_str ' (' info.contacts ')']; end
        end
        title(title_str); ylabel('Amplitude'); grid on;
        if i_ch_plot == length(channel_cols_to_plot), xlabel('Time (s from session start)'); end
    end
    sgtitle(fig_vis, sprintf('Broadband Signals: %s - %s', patient_folder, session_name), 'Interpreter', 'none');
    [~, device_name_short] = fileparts(device_path_sess);
    vis_filename_png = fullfile(vis_dir, sprintf('Vis_%s_%s_%s_%s.png', patient_folder, session_name, device_name_short, 'FirstXMin'));
    try saveas(fig_vis, vis_filename_png); fprintf('  Visualization saved: %s\n', vis_filename_png);
    catch save_err, fprintf('Error saving visualization: %s\n',save_err.message); end
    close(fig_vis);
end
