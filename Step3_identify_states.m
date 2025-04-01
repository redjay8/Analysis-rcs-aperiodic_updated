%% To do:
% Confirm how to identify medication states from eventLogTable
% Attempt differential 1/f biomarker analysis

%% Step3_identify_states.m
% Identify medication and stimulation states from RC+S data
clear;
close all;
clc;

%% 1. Setup paths
base_path = '/../data/';
addpath(genpath('../code/Analysis-rcs-data'));

%% 2. Load combined data
pathName = fullfile(base_path, 'Patient1', 'Device', 'Session');
combined_file = fullfile(pathName, 'combinedDataTable.mat');
fprintf('Loading combined data from: %s\n', combined_file);
load(combined_file);

%% 3. Identify medication states from eventLogTable [TENTATIVE needs testing]
fprintf('Identifying medication states from event log...\n');

%% 3. Identify medication states from eventLogTable
fprintf('Identifying medication states from event log...\n');

% Initialize array to track medication state (0 = OFF, 1 = ON)
med_state = zeros(size(combinedDataTable, 1), 1);

if ~isempty(eventLogTable)
    % Find medication-related events
    med_events_idx = find(contains(eventLogTable.EventType, 'Medication') | ...
                     contains(eventLogTable.EventSubType, 'Medication'));
    
    if ~isempty(med_events_idx)
        fprintf('Found %d medication events in event log\n', length(med_events_idx));
        
        % Process each medication event - assume medication ON for 2 hours after each event
        for i = 1:length(med_events_idx)
            event_idx = med_events_idx(i);
            event_time = eventLogTable.HostUnixTime(event_idx);
            
            % Standard duration of medication effect (2 hours)
            med_duration = 2 * 60 * 60; % 2 hours in seconds
            
            % Assume medication is ON for the duration after each event
            med_state_value = 1; 
            
            % Mark the medication state
            idx = combinedDataTable.DerivedTime >= event_time & ...
                  combinedDataTable.DerivedTime < (event_time + med_duration);
            med_state(idx) = med_state_value;
            
            fprintf('Event %d at %s, duration: %.1f hours, state: %d\n', ...
                   i, datestr(datetime(event_time, 'ConvertFrom', 'posixtime')), ...
                   med_duration/3600, med_state_value);
        end
    else
        fprintf('No medication events found in event log. All data will be treated as OFF medication.\n');
    end
else
    fprintf('Event log is empty. All data will be treated as OFF medication.\n');
end

%% 4. Identify stimulation states [TENTATIVE needs testing]
fprintf('Identifying stimulation states...\n');

% Initialize array to track stimulation state (0 = OFF/low, 1 = ON/high)
stim_state = zeros(size(combinedDataTable, 1), 1);

% Check if using adaptive stimulation
using_adaptive = false;
if any(strcmp(combinedDataTable.Properties.VariableNames, 'Adaptive_CurrentAdaptiveState'))
    % Check if adaptive data is meaningful (i.e., not all NaN)
    if sum(~isnan(combinedDataTable.Adaptive_CurrentAdaptiveState)) > 0
        using_adaptive = true;
        fprintf('Using adaptive stimulation data to determine stimulation states\n');
    end
end

if using_adaptive
    % Get information about adaptive states and their configurations
    high_stim_states = [];
    
    if exist('AdaptiveStimSettings', 'var') && ~isempty(AdaptiveStimSettings)
        % Find the latest adaptive settings entry
        [~, latest_idx] = max(AdaptiveStimSettings.HostUnixTime);
        
        % Get state amplitudes
        state_fields = fieldnames(AdaptiveStimSettings.states);
        
        % Identify high stimulation states (those with amplitudes above threshold)
        stim_threshold = 1.0; % Adjust based on your protocol
        
        for i = 1:length(state_fields)
            state_field = state_fields{i};
            if contains(state_field, 'state')
                state_num = str2double(regexprep(state_field, 'state', ''));
                state_amp = AdaptiveStimSettings.states.(state_field).ampInMilliamps(latest_idx);
                
                if ~isnan(state_amp) && state_amp > stim_threshold
                    high_stim_states = [high_stim_states, state_num];
                end
            end
        end
        
        fprintf('Identified high stimulation states: %s\n', mat2str(high_stim_states));
    else
        % If no settings available, use a heuristic approach
        % Typically higher numbered states have higher stimulation
        high_stim_states = 5:8;
        fprintf('No adaptive settings found. Using default high stimulation states: %s\n', mat2str(high_stim_states));
    end
    
    % Mark times when adaptive state is in high stim
    for i = 1:height(combinedDataTable)
        if ~isnan(combinedDataTable.Adaptive_CurrentAdaptiveState(i)) && ...
           ismember(combinedDataTable.Adaptive_CurrentAdaptiveState(i), high_stim_states)
            stim_state(i) = 1;
        end
    end
    
else
    % Use stimSettingsOut for conventional stimulation
    if ~isempty(stimSettingsOut)
        % Define threshold for "high" stimulation amplitude
        high_stim_threshold = 2.0; % [TENTATIVE]
        
        for i = 1:height(stimSettingsOut)
            % Determine time range for this setting
            if i < height(stimSettingsOut)
                idx = combinedDataTable.DerivedTime >= stimSettingsOut.HostUnixTime(i) & ...
                      combinedDataTable.DerivedTime < stimSettingsOut.HostUnixTime(i+1);
            else
                idx = combinedDataTable.DerivedTime >= stimSettingsOut.HostUnixTime(i);
            end
            
            % Check if stimulation is ON
            if strcmp(stimSettingsOut.therapyStatus{i}, 'On')
                % Get current active group settings
                active_group = stimSettingsOut.activeGroup{i};
                
                switch active_group
                    case 'GroupA'
                        amp = stimSettingsOut.GroupA.ampInMilliamps(i);
                    case 'GroupB'
                        amp = stimSettingsOut.GroupB.ampInMilliamps(i);
                    case 'GroupC'
                        amp = stimSettingsOut.GroupC.ampInMilliamps(i);
                    case 'GroupD'
                        amp = stimSettingsOut.GroupD.ampInMilliamps(i);
                    otherwise
                        amp = 0;
                end
                
                % Mark as high stimulation if amplitude > threshold
                if amp > high_stim_threshold
                    stim_state(idx) = 1;
                end
            end
        end
        
        fprintf('Identified stimulation states using stimSettingsOut\n');
    else
        fprintf('No stimulation data available. All data will be treated as low stimulation.\n');
    end
end

%% 5. Create condition indices
off_med_off_stim_idx = med_state == 0 & stim_state == 0;
off_med_on_stim_idx = med_state == 0 & stim_state == 1;
on_med_off_stim_idx = med_state == 1 & stim_state == 0;
on_med_on_stim_idx = med_state == 1 & stim_state == 1;

fprintf('OFF medication, OFF/low stimulation: %d samples (%.2f%%)\n', 
    sum(off_med_off_stim_idx), 100*sum(off_med_off_stim_idx)/length(med_state));
fprintf('OFF medication, ON/high stimulation: %d samples (%.2f%%)\n', 
    sum(off_med_on_stim_idx), 100*sum(off_med_on_stim_idx)/length(med_state));
fprintf('ON medication, OFF/low stimulation: %d samples (%.2f%%)\n', 
    sum(on_med_off_stim_idx), 100*sum(on_med_off_stim_idx)/length(med_state));
fprintf('ON medication, ON/high stimulation: %d samples (%.2f%%)\n', 
    sum(on_med_on_stim_idx), 100*sum(on_med_on_stim_idx)/length(med_state));

%% 6. Save state information
state_file = fullfile(pathName, 'condition_states.mat');
save(state_file, 'med_state', 'stim_state', 'off_med_off_stim_idx', 'off_med_on_stim_idx', ...
    'on_med_off_stim_idx', 'on_med_on_stim_idx');
fprintf('State information saved to: %s\n', state_file);

fprintf('Step 3 processing complete.\n');