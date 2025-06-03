# Analysis-rcs-aperiodic  

## Step1.m:  
Implements ProcessRCS

## Step2.m:
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

## Step3.ipynb:

#### Data Ingestion & Clinical State Integration: 
Loading pre-calculated PSDs and associated metadata for a specific subject and hemisphere. Then integrate clinical state information (e.g., mobile, immobile, sleep, dyskinesia levels based on PKG scores) using a 2-minute windowing approach to assign a dominant clinical state to each neural data segment.

#### Spectral Feature Extraction with FOOOF:
Core of the analysis involves using the FOOOF (Fitting Oscillations & One Over F) algorithm. For each individual PSD segment, FOOOF is applied across multiple predefined frequency ranges (e.g., 10-40Hz, 30-90Hz, 10-90Hz) using two aperiodic settings: 'fixed' (1/f) and 'knee' (1/f with a knee parameter). This separates the aperiodic (1/f-like) background activity from true periodic oscillations (peaks).
Key aperiodic parameters (offset, exponent, and knee if applicable), along with goodness-of-fit metrics (R-squared, error) and peak parameters, are extracted for each fit.
The script identifies "oscillatory humps"—broad regions of power rising above the fitted aperiodic component—and characterizes their frequency range, width, and maximum power.

#### Targeted Band Power Analysis:
To analyze specific frequency bands like beta (13-30Hz) and gamma (60-90Hz), the script first determines a channel-specific "dominant" peak frequency within these bands. This is achieved by:
            Temporarily fitting a FOOOF model to each segment to estimate its aperiodic component.
            Subtracting this aperiodic fit from the log-power spectrum to create a "flattened" spectrum.
            Identifying significant peaks in the flattened spectrum within the beta and gamma bands, based on a threshold relative to a baseline power region (e.g., 10-12Hz for beta, 55-59Hz for gamma).
            The mode of these significant peak frequencies across all segments for a given channel is deemed its dominant beta/gamma frequency.
        Finally, the power at these channel-specific dominant beta and gamma frequencies is extracted from the original (non-flattened) log-power spectra of each segment. This approach aims to capture the power of the most prominent, consistent oscillation within the band, rather than an average power that might be diluted by aperiodic shifts or broader activity.

#### Averaging & Visualization:
The script calculates and visualizes average PSDs for each recording channel, grouped by the derived clinical states. FOOOF is then run on these averaged PSDs to examine state-dependent aperiodic characteristics.
It also generates plots for individual FOOOF fits and distributions of oscillatory hump widths.

#### Master Data Compilation: 
All extracted features—including original segment metadata, clinical states, LEDD, dominant beta/gamma peak powers, detailed parameters from both 'fixed' and 'knee' FOOOF fits for each frequency band, and a "best model" selection based on R-squared—are compiled into a comprehensive master CSV file for subsequent statistical analysis.

This systematic approach allows for a detailed characterization of how both aperiodic and periodic neural activity relate to varying clinical states and medication levels. The use of FOOOF is critical for disentangling these two components of the neural power spectrum, which builds the foundation of investigating whether aperiodics offer orthogonal value over beta and gamma power.


