# Analysis-rcs-aperiodic  

## Descriptions by script:

### Step1.m:  
Implements ProcessRCS

### Step2.m:
Performs center-aligned temporal segmentation of neural time series data with contralateral PKG (Parkinson's Kinetigraph) measurements for spectral analysis. Key processing steps include:
Data Integration:

Loads preprocessed RC+S neural recordings from Step 1 master data files containing multichannel time series
Retrieves contralateral PKG accelerometry scores (BK, DK, tremor metrics) to capture motor symptom severity
Implements UTC-aware timestamp alignment between neural and behavioral data streams

Signal Preprocessing:

Applies 4th-order Butterworth high-pass filter (1 Hz cutoff) to remove DC drift and slow physiological artifacts
Preserves pathological low-frequency oscillations (4-30 Hz) critical for Parkinsonian biomarkers
Segments only continuous data chunks >5 seconds to ensure filter stability

Temporal Alignment Strategy:

Centers 120-second neural windows around 30-second interpolated PKG timepoints
This 4:1 ratio balances spectral resolution (0.5 Hz with 2s Welch windows) with behavioral state stability
Center alignment minimizes temporal mismatch between neural activity and motor symptoms

Spectral Analysis:

Computes PSDs using Welch's method (2s Hanning windows, 50% overlap) for robust spectral estimation
Parallel processing across channels maximizes computational efficiency
Outputs include full spectral data for FOOOF decomposition and band power analyses



### Step4b_visualize_fooof_fits.m
Now able to visualize individual FOOOF fits for selected time windows from LFP data.
% NOTE: only validated on demo data as of now.

### Step4_calculate_rolling_aperiodic.m  
Now adopting three methods accounting for mismatch in neural data and PKG sampling frequencies  
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
