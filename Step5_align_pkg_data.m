%% Step5_align_pkg_data.m
% Import PKG data from CSV format and align with aperiodic exponents
% Performs sensitivity testing across all three methods
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));

%% 2. Load aperiodic data from all methods and condition states
patient = 'Patient1';  % RC+S recording patient name
pkg_patient = 'Pat1';  % Corresponding PKG patient name (if different)
side = 'L';
brain_regions = {'+2-0', '+3-1', '+9-8', '+11-10'}; 
brain_region_name = {'STN1', 'STN1', 'postcentral gyrus', 'precentral gyrus'}; 
chan = 2;

% Path name
pathName = fullfile(base_path, patient, 'Device', 'Session');

% Load condition states (medication and stimulation states)
state_file = fullfile(pathName, 'condition_states.mat');
fprintf('Loading condition states from: %s\n', state_file);
load(state_file);

% Load aperiodic exponent time series results from all three methods
all_methods_file = fullfile(pathName, sprintf('aperiodic_timeseries_all_%s_%s.mat', brain_regions{chan}, side));
if exist(all_methods_file, 'file')
    fprintf('Loading all aperiodic methods from: %s\n', all_methods_file);
    load(all_methods_file);
    all_methods_loaded = true;
else
    fprintf('Combined methods file not found. Looking for individual files...\n');
    all_methods_loaded = false;
    
    % Try to load Method 1
    m1_file = fullfile(pathName, sprintf('aperiodic_timeseries_%s_%s.mat', brain_regions{chan}, side));
    if exist(m1_file, 'file')
        fprintf('Loading Method 1 from: %s\n', m1_file);
        temp = load(m1_file);
        aperiodic_results_m1 = temp.aperiodic_results;
        metadata = temp.metadata;
        m1_loaded = true;
    else
        error('No aperiodic results files found. Run Step 4 first.');
    end
    
    % Try to load Method 2
    m2_file = fullfile(pathName, sprintf('aperiodic_timeseries_m2_%s_%s.mat', brain_regions{chan}, side));
    if exist(m2_file, 'file')
        fprintf('Loading Method 2 from: %s\n', m2_file);
        temp = load(m2_file);
        aperiodic_results_m2 = temp.aperiodic_results;
        m2_loaded = true;
    else
        m2_loaded = false;
        fprintf('Method 2 file not found.\n');
    end
    
    % Try to load Method 3
    m3_file = fullfile(pathName, sprintf('aperiodic_timeseries_m3_%s_%s.mat', brain_regions{chan}, side));
    if exist(m3_file, 'file')
        fprintf('Loading Method 3 from: %s\n', m3_file);
        temp = load(m3_file);
        aperiodic_results_m3 = temp.aperiodic_results;
        m3_loaded = true;
    else
        m3_loaded = false;
        fprintf('Method 3 file not found.\n');
    end
end

%% 3. Load PKG data from CSV file(s)
fprintf('Loading PKG data...\n');

% PKG data is stored in CSV format with data every 2 minutes (120 seconds)
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

fprintf('Loaded %d PKG data points across %d files\n', height(pkg_data), length(pkg_files));

%% 4. Process PKG data
fprintf('Processing PKG data...\n');

% Convert date-time to POSIX time for easier alignment
pkg_data.datetime = datetime(pkg_data.Date_Time);
pkg_data.time = posixtime(pkg_data.datetime);

% Remove periods when PKG was off wrist
if ismember('Off_Wrist', pkg_data.Properties.VariableNames)
    pkg_data = pkg_data(~pkg_data.Off_Wrist, :);
    fprintf('Removed off-wrist periods. %d data points remaining.\n', height(pkg_data));
end

% Sort by time for proper alignment
pkg_data = sortrows(pkg_data, 'time');

% Display statistics about PKG symptom measures
fprintf('PKG data summary statistics:\n');
if ismember('BKS', pkg_data.Properties.VariableNames)
    fprintf('BKS: Mean=%.2f, SD=%.2f, Min=%.2f, Max=%.2f\n', ...
        mean(pkg_data.BKS), std(pkg_data.BKS), min(pkg_data.BKS), max(pkg_data.BKS));
end
if ismember('DKS', pkg_data.Properties.VariableNames)
    fprintf('DKS: Mean=%.2f, SD=%.2f, Min=%.2f, Max=%.2f\n', ...
        mean(pkg_data.DKS), std(pkg_data.DKS), min(pkg_data.DKS), max(pkg_data.DKS));
end
if ismember('TS', pkg_data.Properties.VariableNames)
    fprintf('TS: Mean=%.2f, SD=%.2f, Min=%.2f, Max=%.2f\n', ...
        mean(pkg_data.TS), std(pkg_data.TS), min(pkg_data.TS), max(pkg_data.TS));
end

%% 5. Align PKG data with aperiodic exponents for all methods
fprintf('Aligning PKG data with aperiodic exponents from all methods...\n');

% --- Method 1: Original method (5s steps) ---
if exist('aperiodic_results_m1', 'var')
    fprintf('Aligning Method 1 (5s steps)...\n');
    aligned_data_m1 = align_pkg_with_aperiodic(aperiodic_results_m1, pkg_data);
    
    % Calculate alignment statistics
    pkg_alignment_rate_m1 = 100 * sum(~isnan(aligned_data_m1.pkg_BKS)) / height(aligned_data_m1);
    fprintf('Method 1: %.2f%% of aperiodic data points have corresponding PKG data\n', pkg_alignment_rate_m1);
end

% --- Method 2: 60s windows ending at PKG times ---
if exist('aperiodic_results_m2', 'var')
    fprintf('Aligning Method 2 (aligned to PKG times)...\n');
    aligned_data_m2 = align_pkg_with_aperiodic(aperiodic_results_m2, pkg_data);
    
    % Calculate alignment statistics
    pkg_alignment_rate_m2 = 100 * sum(~isnan(aligned_data_m2.pkg_BKS)) / height(aligned_data_m2);
    fprintf('Method 2: %.2f%% of aperiodic data points have corresponding PKG data\n', pkg_alignment_rate_m2);
end

% --- Method 3: 120s window averaging ---
if exist('aperiodic_results_m3', 'var')
    fprintf('Aligning Method 3 (120s window averaging)...\n');
    aligned_data_m3 = align_pkg_with_aperiodic(aperiodic_results_m3, pkg_data);
    
    % Calculate alignment statistics
    pkg_alignment_rate_m3 = 100 * sum(~isnan(aligned_data_m3.pkg_BKS)) / height(aligned_data_m3);
    fprintf('Method 3: %.2f%% of aperiodic data points have corresponding PKG data\n', pkg_alignment_rate_m3);
end

%% 6. Save aligned data for all methods
fprintf('Saving aligned data for all methods...\n');
aligned_data_file = fullfile(pathName, sprintf('aligned_pkg_aperiodic_all_%s_%s.mat', brain_regions{chan}, side));

% Make sure we save the metadata
if ~exist('metadata', 'var')
    metadata = struct();
    metadata.patient = patient;
    metadata.side = side;
    metadata.brain_region = brain_regions{chan};
    metadata.brain_region_name = brain_region_name{chan};
end

save(aligned_data_file, 'aligned_data_m1', 'aligned_data_m2', 'aligned_data_m3', 'metadata', 'pkg_data');
fprintf('Aligned data for all methods saved to: %s\n', aligned_data_file);

%% 7. Perform sensitivity analysis across methods
fprintf('Performing sensitivity analysis across methods...\n');

% Create a figure to compare correlation strengths across methods
figure('Position', [100 100 1500 800]);

% Define symptom measures to analyze
measures = {'BKS', 'DKS', 'TS'};
method_names = {'Method 1: 5s steps', 'Method 2: 60s ending at PKG', 'Method 3: 120s avg before PKG'};
method_colors = {'b', 'r', 'g'};

% Prepare a table for storing all correlation results
corr_results = table();
row_idx = 1;

% Iterate through each symptom measure
for m = 1:length(measures)
    measure = measures{m};
    subplot(1, length(measures), m);
    hold on;
    
    % Store correlation values for comparison
    corr_values = [];
    p_values = [];
    method_indices = [];
    
    % Analyze each method if available
    methods = {'aligned_data_m1', 'aligned_data_m2', 'aligned_data_m3'};
    
    for a = 1:length(methods)
        method = methods{a};
        if exist(method, 'var')
            % Get data for this method
            data = eval(method);
            
            % Check if PKG measure exists
            pkg_col = ['pkg_' measure];
            if ~ismember(pkg_col, data.Properties.VariableNames)
                fprintf('PKG measure %s not found in %s\n', measure, method_names{a});
                continue;
            end
            
            % Compute correlation
            valid_idx = ~isnan(data.(pkg_col)) & ~isnan(data.exponent);
            
            if sum(valid_idx) > 5
                [r, p] = corrcoef(data.(pkg_col)(valid_idx), data.exponent(valid_idx));
                corr_val = r(1,2);
                p_val = p(1,2);
                
                % Store for comparison
                corr_values = [corr_values; corr_val];
                p_values = [p_values; p_val];
                method_indices = [method_indices; a];
                
                % Add to correlation results table
                corr_results.Method(row_idx) = {method_names{a}};
                corr_results.Measure(row_idx) = {measure};
                corr_results.Correlation(row_idx) = corr_val;
                corr_results.PValue(row_idx) = p_val;
                corr_results.N(row_idx) = sum(valid_idx);
                row_idx = row_idx + 1;
                
                % Plot scatter with regression line
                scatter(data.(pkg_col)(valid_idx), data.exponent(valid_idx), 30, method_colors{a}, 'filled', 'MarkerFaceAlpha', 0.3);
                
                % Add regression line
                p_fit = polyfit(data.(pkg_col)(valid_idx), data.exponent(valid_idx), 1);
                x_range = linspace(min(data.(pkg_col)(valid_idx)), max(data.(pkg_col)(valid_idx)), 100);
                y_fit = polyval(p_fit, x_range);
                plot(x_range, y_fit, method_colors{a}, 'LineWidth', 2);
                
                % Add correlation value text with position adjustment to avoid overlap
                text_x = max(data.(pkg_col)(valid_idx))*0.7;
                text_y = max(data.exponent(valid_idx))*(0.9 - 0.1*a);
                text_str = sprintf('%s: r=%.3f, p=%.3f, n=%d', method_names{a}, corr_val, p_val, sum(valid_idx));
                text(text_x, text_y, text_str, 'Color', method_colors{a}, 'FontWeight', 'bold');
            else
                fprintf('Not enough valid data points for %s in %s\n', measure, method_names{a});
            end
        end
    end
    
    % Add plot details
    title(sprintf('%s vs. Aperiodic Exponent - Method Comparison', measure));
    xlabel(sprintf('%s Score', measure));
    ylabel('Aperiodic Exponent');
    grid on;
    
    % Only add legend if we have multiple methods
    if length(method_indices) > 1
        methods_plotted = method_names(method_indices);
        legend(methods_plotted, 'Location', 'best');
    end
end

% Save correlation results table
corr_results_file = fullfile(pathName, sprintf('correlation_results_%s_%s.mat', brain_regions{chan}, side));
save(corr_results_file, 'corr_results');

% Display table
fprintf('\nCorrelation Results for all methods:\n');
disp(corr_results);

% Save figure
saveas(gcf, fullfile(pathName, sprintf('sensitivity_analysis_%s_%s.fig', brain_regions{chan}, side)));
saveas(gcf, fullfile(pathName, sprintf('sensitivity_analysis_%s_%s.png', brain_regions{chan}, side)));

%% 8. Medication and Stimulation State Analysis across methods
fprintf('Analyzing effects of medication and stimulation states across methods...\n');

figure('Position', [100 100 1500 800]);

% Define states to analyze
states = {'med_state', 'stim_state'};
state_labels = {'Medication State', 'Stimulation State'};
state_values = {[0, 1], [0, 1]};
state_value_labels = {{'OFF', 'ON'}, {'OFF/Low', 'ON/High'}};

% State comparison results table
state_results = table();
state_row = 1;

% Iterate through methods
for a = 1:length(methods)
    method = methods{a};
    if exist(method, 'var')
        data = eval(method);
        
        for s = 1:length(states)
            state = states{s};
            subplot(length(methods), length(states), (a-1)*length(states) + s);
            
            % Group by state values
            state_groups = {};
            exp_values = [];
            
            for v = 1:length(state_values{s})
                state_val = state_values{s}(v);
                value_idx = data.(state) == state_val & ~isnan(data.exponent);
                
                if sum(value_idx) > 0
                    state_groups = [state_groups; repmat({state_value_labels{s}{v}}, sum(value_idx), 1)];
                    exp_values = [exp_values; data.exponent(value_idx)];
                end
            end
            
            % Create box plot if data is available
            if ~isempty(state_groups)
                boxplot(exp_values, state_groups);
                title(sprintf('%s by %s - %s', 'Aperiodic Exponent', state_labels{s}, method_names{a}));
                ylabel('Aperiodic Exponent');
                grid on;
                
                % Perform t-test if there are two groups with sufficient data
                unique_groups = unique(state_groups);
                if length(unique_groups) == 2
                    group1_idx = strcmp(state_groups, unique_groups{1});
                    group2_idx = strcmp(state_groups, unique_groups{2});
                    
                    if sum(group1_idx) > 5 && sum(group2_idx) > 5
                        [h, p, ~, stats] = ttest2(exp_values(group1_idx), exp_values(group2_idx));
                        
                        % Store in state results table
                        state_results.Method(state_row) = {method_names{a}};
                        state_results.State(state_row) = {state_labels{s}};
                        state_results.Group1(state_row) = unique_groups(1);
                        state_results.Group2(state_row) = unique_groups(2);
                        state_results.Mean1(state_row) = mean(exp_values(group1_idx));
                        state_results.Mean2(state_row) = mean(exp_values(group2_idx));
                        state_results.P_Value(state_row) = p;
                        state_results.T_Stat(state_row) = stats.tstat;
                        state_results.DF(state_row) = stats.df;
                        state_row = state_row + 1;
                        
                        fprintf('%s - %s: t(%d)=%.3f, p=%.4f, %s vs %s\n', ...
                            method_names{a}, state_labels{s}, stats.df, stats.tstat, p, ...
                            unique_groups{1}, unique_groups{2});
                        
                        % Add p-value to plot
                        text(1.5, max(exp_values)*0.9, sprintf('p = %.4f', p), ...
                            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                    end
                end
            else
                text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center');
                axis off;
            end
        end
    end
end

% Save state results table
state_results_file = fullfile(pathName, sprintf('state_results_%s_%s.mat', brain_regions{chan}, side));
save(state_results_file, 'state_results');

% Display table
fprintf('\nState Comparison Results for all methods:\n');
disp(state_results);

% Save figure
saveas(gcf, fullfile(pathName, sprintf('state_analysis_%s_%s.fig', brain_regions{chan}, side)));
saveas(gcf, fullfile(pathName, sprintf('state_analysis_%s_%s.png', brain_regions{chan}, side)));

%% 9. Time-series visualization of best method (among the three) with PKG data
fprintf('Creating time-series visualization of aperiodic exponents with PKG data...\n');

% Determine which method has the strongest correlation
best_method_idx = 1; % Default to Method 1
if exist('corr_results', 'var') && height(corr_results) > 0
    % Find method with strongest average correlation across measures
    method_avg_corr = zeros(length(methods), 1);
    method_exists = false(length(methods), 1);
    
    for a = 1:length(methods)
        method = method_names{a};
        method_rows = strcmp(corr_results.Method, method);
        
        if any(method_rows)
            method_avg_corr(a) = mean(abs(corr_results.Correlation(method_rows)));
            method_exists(a) = true;
        end
    end
    
    % Select best method from those that exist
    [~, best_among_existing] = max(method_avg_corr(method_exists));
    existing_methods = find(method_exists);
    best_method_idx = existing_methods(best_among_existing);
end

% Get best method data
best_method = methods{best_method_idx};
if exist(best_method, 'var')
    best_data = eval(best_method);
    
    figure('Position', [100 100 1200 600]);
    
    % Plot aperiodic exponent time series
    ax1 = subplot(2,1,1);
    plot(best_data.datetime, best_data.exponent, 'LineWidth', 1.5);
    title(sprintf('Aperiodic Exponent Time Series (%s)', method_names{best_method_idx}));
    ylabel('Exponent');
    grid on;
    
    % Highlight by medication state
    hold on;
    med_off_idx = best_data.med_state == 0 & ~isnan(best_data.exponent);
    med_on_idx = best_data.med_state == 1 & ~isnan(best_data.exponent);
    scatter(best_data.datetime(med_off_idx), best_data.exponent(med_off_idx), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    scatter(best_data.datetime(med_on_idx), best_data.exponent(med_on_idx), 30, 'g', 'filled', 'MarkerFaceAlpha', 0.3);
    legend('Exponent', 'OFF Medication', 'ON Medication', 'Location', 'best');
    
    % Plot PKG measures
    ax2 = subplot(2,1,2);
    hold on;
    
    % Plot all available PKG measures
    if ismember('pkg_BKS', best_data.Properties.VariableNames)
        plot(best_data.datetime, best_data.pkg_BKS, 'LineWidth', 1.5, 'Color', 'r');
    end
    if ismember('pkg_DKS', best_data.Properties.VariableNames)
        plot(best_data.datetime, best_data.pkg_DKS, 'LineWidth', 1.5, 'Color', 'g');
    end
    if ismember('pkg_TS', best_data.Properties.VariableNames)
        plot(best_data.datetime, best_data.pkg_TS, 'LineWidth', 1.5, 'Color', 'b');
    end
    
    title('PKG Measures');
    ylabel('Score');
    xlabel('Time');
    grid on;
    
    % Add legend with appropriate measures
    legend_items = {};
    if ismember('pkg_BKS', best_data.Properties.VariableNames)
        legend_items{end+1} = 'BKS';
    end
    if ismember('pkg_DKS', best_data.Properties.VariableNames)
        legend_items{end+1} = 'DKS';
    end
    if ismember('pkg_TS', best_data.Properties.VariableNames)
        legend_items{end+1} = 'TS';
    end
    legend(legend_items, 'Location', 'best');
    
    % Link axes for synchronized zooming
    linkaxes([ax1, ax2], 'x');
    
    % Save figure
    saveas(gcf, fullfile(pathName, sprintf('timeseries_best_method_%s_%s.fig', brain_regions{chan}, side)));
    saveas(gcf, fullfile(pathName, sprintf('timeseries_best_method_%s_%s.png', brain_regions{chan}, side)));
end

fprintf('Step 5 processing complete with sensitivity analysis.\n');

%% Helper function to align PKG data with aperiodic exponents
function aligned_data = align_pkg_with_aperiodic(aperiodic_results, pkg_data)
    % Convert aperiodic times to datetime for easier processing
    aperiodic_datetime = datetime(aperiodic_results.time, 'ConvertFrom', 'posixtime');

    % Create aligned table
    aligned_data = table();
    aligned_data.time = aperiodic_results.time;
    aligned_data.datetime = aperiodic_datetime;
    aligned_data.exponent = aperiodic_results.exponent;
    aligned_data.offset = aperiodic_results.offset;
    aligned_data.r_squared = aperiodic_results.r_squared;
    aligned_data.med_state = aperiodic_results.med_state;
    aligned_data.stim_state = aperiodic_results.stim_state;

    % Initialize PKG fields in aligned table with all potential measures
    if ismember('BKS', pkg_data.Properties.VariableNames)
        aligned_data.pkg_BKS = NaN(height(aligned_data), 1);
    end
    if ismember('DKS', pkg_data.Properties.VariableNames)
        aligned_data.pkg_DKS = NaN(height(aligned_data), 1);
    end
    if ismember('TS', pkg_data.Properties.VariableNames)
        aligned_data.pkg_TS = NaN(height(aligned_data), 1);
    end

    % For each aperiodic data point, find the closest FUTURE PKG data point
    for i = 1:height(aligned_data)
        % Get current time
        current_time = aligned_data.time(i);
        
        % Find times that are AFTER the current time (future PKG readings)
        future_indices = find(pkg_data.time > current_time);
        
        if ~isempty(future_indices)
            % Find the closest future PKG data point
            [~, min_idx] = min(pkg_data.time(future_indices) - current_time);
            closest_idx = future_indices(min_idx);
            
            % Calculate time difference in seconds
            time_diff = pkg_data.time(closest_idx) - current_time;
            
            % Only use if within 5 minutes (300 seconds)
            if time_diff <= 300
                % Use consistent checks for all PKG measures
                if ismember('BKS', pkg_data.Properties.VariableNames) && ismember('pkg_BKS', aligned_data.Properties.VariableNames)
                    aligned_data.pkg_BKS(i) = pkg_data.BKS(closest_idx);
                end
                
                if ismember('DKS', pkg_data.Properties.VariableNames) && ismember('pkg_DKS', aligned_data.Properties.VariableNames)
                    aligned_data.pkg_DKS(i) = pkg_data.DKS(closest_idx);
                end
                
                if ismember('TS', pkg_data.Properties.VariableNames) && ismember('pkg_TS', aligned_data.Properties.VariableNames)
                    aligned_data.pkg_TS(i) = pkg_data.TS(closest_idx);
                end
            end
        end
    end
end