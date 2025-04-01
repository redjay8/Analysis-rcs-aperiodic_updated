function [ap_params] = extract_aperiodic_parameters(freqs, power, freq_range, settings)
    % EXTRACT_APERIODIC_PARAMETERS Extract aperiodic parameters using a simplified FOOOF approach
    %
    % Inputs:
    %   freqs - Frequency vector
    %   power - Power spectrum
    %   freq_range - Frequency range for analysis [min_freq max_freq]
    %   settings - FOOOF settings structure
    %
    % Outputs:
    %   ap_params - Structure containing aperiodic parameters
    
    % Initialize output structure
    ap_params = struct();
    
    % Find indices for specified frequency range
    freq_idx = freqs >= freq_range(1) & freqs <= freq_range(2);
    ap_params.freq_idx = freq_idx;
    
    % Extract relevant frequencies and power
    freqs_fit = freqs(freq_idx);
    power_fit = power(freq_idx);
    
    % Take log of power and frequency for linear fitting
    log_freqs = log10(freqs_fit);
    log_power = log10(power_fit);
    
    % Initial parameters for aperiodic fit: [offset, exponent]
    % Estimate using linear regression in log-log space
    p = polyfit(log_freqs, log_power, 1);
    init_params = [mean(log_power) - p(1)*mean(log_freqs), -p(1)];
    
    % Define aperiodic function (1/f^χ in log-log space)
    ap_func = @(p, f) p(1) - p(2)*log10(f);
    
    % Fit aperiodic component
    options = optimset('Display', 'off');
    [ap_params_fit, ~, ~, ~, ~, ~, ~] = lsqcurvefit(ap_func, init_params, freqs_fit, log_power, [], [], options);
    
    % Extract parameters
    ap_params.offset = ap_params_fit(1);
    ap_params.exponent = ap_params_fit(2);
    
    % Generate aperiodic fit
    ap_params.aperiodic_fit = 10.^ap_func(ap_params_fit, freqs_fit);
    
    % Calculate goodness of fit
    ss_total = sum((log_power - mean(log_power)).^2);
    ss_residual = sum((log_power - ap_func(ap_params_fit, freqs_fit)).^2);
    ap_params.r_squared = 1 - (ss_residual / ss_total);
end