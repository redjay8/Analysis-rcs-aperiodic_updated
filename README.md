# Analysis-rcs-aperiodic  

## Descriptions by script:

### Step1.m:  
Implements ProcessRCS

### Step2.m:
Performs center-aligned temporal segmentation of neural time series data with contralateral PKG (Parkinson's Kinetigraph) measurements for spectral analysis. Key processing steps include:

#### Data Integration:

Loads preprocessed RC+S neural recordings from Step 1 master data files containing multichannel time series
Retrieves contralateral PKG accelerometry scores (BK, DK, tremor metrics) to capture motor symptom severity
Implements UTC-aware timestamp alignment between neural and behavioral data streams

#### Signal Preprocessing:

Applies 4th-order Butterworth high-pass filter (1 Hz cutoff) to remove DC drift and slow physiological artifacts
Preserves pathological low-frequency oscillations (4-30 Hz) critical for Parkinsonian biomarkers
Segments only continuous data chunks >5 seconds to ensure filter stability

#### Temporal Alignment Strategy:

Centers 120-second neural windows around 30-second interpolated PKG timepoints
This 4:1 ratio balances spectral resolution (0.5 Hz with 2s Welch windows) with 2-min updating PKG data
Center alignment minimizes temporal mismatch between neural activity and motor symptoms

#### Spectral Analysis:

Computes PSDs using Welch's method (2s Hanning windows, 50% overlap) for robust spectral estimation
Parallel processing across channels for computational efficiency
Outputs include full spectral data for FOOOF decomposition and band power analyses
