
function Batch_Step2_Runx()

    clear; close all; clc;
    % === Config you can tweak if needed ===
    root_dir      = '/home/jackson/Downloads/Folder';                      % <-- your actual root
    toolbox_path  = '/home/jackson/Analysis-rcs-data';
    % If Step2_github.m is not on path, hardcode it here:
    fallback_step2 = '/home/jackson/Analysis-rcs-aperiodic-main/Step2_github.m';

    % Ensure toolbox is available
    if isfolder(toolbox_path)
        addpath(genpath(toolbox_path));
    else
        warning('Toolbox path not found: %s', toolbox_path);
    end

    % Find Step2 script
    step2_script = which('Step2_github.m');
    if isempty(step2_script)
        step2_script = fallback_step2;
        if ~isfile(step2_script)
            error('Could not locate Step2_github.m. Put it on the path or update fallback_step2.');
        end
    end

    patients = {'RCS05'};
    sides    = {'Left'};              % long-form for filenames
    sideLR   = containers.Map({'Left','Right'},{'L','R'});

    fprintf('\n=== Batch Step 2 runs (sequential) ===\nRoot: %s\nScript: %s\n\n', root_dir, step2_script);

    % Keep memory footprint small: use threads pool for parfor (optional)
    pool = gcp('nocreate');
    if ~isempty(pool), delete(pool); end
    try
        parpool('threads');
    catch MEp
        fprintf('Note: could not start thread pool (%s). Proceeding without it.\n', MEp.message);
    end

    for ip = 1:numel(patients)
        patient = patients{ip};
        pat_dir = fullfile(root_dir, patient);

        if ~isfolder(pat_dir)
            fprintf('Skipping %s (folder not found): %s\n', patient, pat_dir);
            continue;
        end

        for iside = 1:numel(sides)
            side = sides{iside};      % 'Left' or 'Right'
            LR   = sideLR(side);      % 'L' or 'R'

            % Expected Step 1 file pattern per your layout:
            mat_name = sprintf('Step1_MasterData_%s%s_%s_AllSessions.mat', patient, LR, side);
            mat_path = fullfile(pat_dir, mat_name);

            if ~isfile(mat_path)
                fprintf('⏭️  %s %s: Step 1 file not found (%s) -> skipping\n', patient, side, mat_name);
                continue;
            end

            % Build batch config to pass into Step2_github.m
            BATCH_MODE = true; %#ok<NASGU>
            outdir     = fullfile(pat_dir, 'Step2newx', sprintf('%s%s_%s', patient, LR, side));
            if ~isfolder(outdir), mkdir(outdir); end

            BATCH_CFG = struct( ... %#ok<NASGU>
                'target_patient_folder_step1', sprintf('%s%s', patient, LR), ...
                'target_hemisphere_step1',     side, ...
                'processed_data_folder_step1', pat_dir, ...
                'pkg_base_data_folder',        fullfile(pat_dir, 'PKG'), ...
                'project_base_path',           pat_dir, ...
                'preprocessed_output_folder',  outdir ...
            );

            % Per-run log
            log_file = fullfile(outdir, sprintf('Step2_log_%s_%s.txt', patient, side));
            try
                diary off; diary(log_file);
            catch
                % if diary fails we still proceed
            end

            fprintf('\n==== Running Step 2 for %s %s ====\n', patient, side);
            fprintf('MAT file: %s\n', mat_path);
            fprintf('PKG dir : %s\n', fullfile(pat_dir,'PKG'));
            fprintf('Output  : %s\n', outdir);

            try
                % Run the script in this workspace; it will see BATCH_MODE/BATCH_CFG
                run(step2_script);
                fprintf('✅ Finished: %s %s\n', patient, side);
            catch ME
                fprintf('❌ ERROR in %s %s: %s\n', patient, side, ME.message);
                % Continue to next job
            end

            % Close log
            try, diary off; catch, end

            % Cleanup before the next run (memory-friendly)
            close all force;
            clearvars -except root_dir toolbox_path fallback_step2 step2_script patients sides sideLR ip iside patient pat_dir
            % no 'pack' (unnecessary & noisy on 64-bit)
        end
    end

    fprintf('\n=== Batch complete ===\n');
end
