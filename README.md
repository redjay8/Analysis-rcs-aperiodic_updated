# Analysis-rcs-aperiodic  

### %% To do:  
Setup differential 1/f biomarker analysis based on identifying medication states from eventLogTable  
Test pipeline on some data  

### %% Changes:  
### Step4_calculate_rolling_aperiodic.m  
Now adopting three methods due to mismatch in neural data and PKG sampling frequencies  
% Method 1 (60s windows with 5s steps):  
% High temporal resolution, but multiple neural measures map to one PKG reading  
% Method 2: Hybrid approach (120s steps, aligned to end at PKG times):  
% Independent measurements with no window overlap  
% Matches the sampling rate of both data streams  
% Method 3: Averaging approach (5s steps, averaged per PKG interval):  
% Leverages all data but reduces dimensionality  
  
### Step5_align_pkg_data.m  
% Framework now complete and ready for testing  
% Import PKG data from CSV format and align with aperiodic exponents  
% Performs sensitivity testing across all three methods  
% Basic summary of exploratory statistics
