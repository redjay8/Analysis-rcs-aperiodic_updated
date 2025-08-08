%% Analyze Contact_x_x Columns in All Session .mat Files
clear; clc;

%% Configuration
folder_path = '/home/jackson/step1_processed_data_multi_session_final_newnaming';
table_var_name = 'combinedDataTable_sess';

mat_files = dir(fullfile(folder_path, '*.mat'));
fprintf('🔍 Found %d MAT files in folder: %s\n', length(mat_files), folder_path);

% Initialize output table
long_summary = table();

%% Loop through session .mat files
for file_idx = 1:length(mat_files)
    mat_file = mat_files(file_idx).name;
    mat_path = fullfile(folder_path, mat_file);
    fprintf('\n📂 [%d/%d] Loading %s\n', file_idx, length(mat_files), mat_file);

    if ~isfile(mat_path)
        warning('❌ File not found: %s', mat_path);
        continue;
    end

    try
        data = load(mat_path);
    catch ME
        warning('⚠️ Error loading %s: %s', mat_path, ME.message);
        continue;
    end

    % Detect and load the session table
    if isfield(data, table_var_name)
        T = data.(table_var_name);
    else
        table_vars = fieldnames(data);
        table_vars = table_vars(structfun(@istable, data));
        if isempty(table_vars)
            warning('❌ No table found in: %s', mat_file);
            continue;
        end
        [~, idx] = max(cellfun(@(x) height(data.(x)), table_vars));
        T = data.(table_vars{idx});
        fprintf('📌 Auto-selected table: %s\n', table_vars{idx});
    end

    total_rows = height(T);
    if total_rows == 0
        warning('⚠️ Table is empty in %s', mat_file);
        continue;
    end

    % Identify Contact_x_x columns
    varnames = T.Properties.VariableNames;
    contact_pattern = '^key\d+_contact_\d+_\d+$';
    contact_columns = varnames(~cellfun('isempty', regexp(varnames, contact_pattern)));
    %contactCols = ~cellfun('isempty', regexp(combinedDataTable_sess.Properties.VariableNames, '^key\d+_contact_\d+_\d+$'));
    if isempty(contact_columns)
        fprintf('⚠️ No Contact_x_x columns found in %s\n', mat_file);
        continue;
    end

    % Compute valid data stats per contact
    for i = 1:length(contact_columns)
        cname = contact_columns{i};
        cdata = T.(cname);
        if ~isnumeric(cdata)
            try
                cdata = str2double(string(cdata));
            catch
                cdata = nan(total_rows,1);
            end
        end
        valid_count = sum(~isnan(cdata));
        pct_valid = 100 * valid_count / total_rows;

        % Append row to summary
        new_row = table(string(mat_file), string(cname), valid_count, total_rows, pct_valid, ...
            'VariableNames', {'Filename', 'Contact', 'ValidCount', 'TotalRows', 'PercentValid'});
        long_summary = [long_summary; new_row];
    end
end

%% Save to CSV
output_csv = fullfile(folder_path, 'contact_validity_by_session.csv');
writetable(long_summary, output_csv);
fprintf('\n📄 Contact validity summary saved to:\n%s\n', output_csv);
