
# SSVEP Target Detection with CCA in MATLAB

MATLAB implementation of an SSVEP-based brain–computer interface (BCI) pipeline for **offline** and **online target detection** using **Canonical Correlation Analysis (CCA)**.

---

# Overview

Steady-State Visual Evoked Potentials (SSVEPs) are brain responses elicited when a subject focuses on a visual stimulus flickering at a specific frequency. In SSVEP-based BCIs, multiple visual targets flicker at different frequencies, and the user selects a target by focusing attention on it.

The EEG signal, therefore, contains frequency components corresponding to the attended stimulus and its harmonics.

This repository implements a CCA-based approach that compares EEG segments with synthetic sinusoidal reference templates and predicts the attended target frequency.

The repository includes both:

- **Offline analysis** (evaluate classification accuracy using recorded trials)
- **Online-style simulation** (process EEG buffers sequentially to emulate real-time detection)

---

# Repository Structure

```
.
├── prepare_offline_data.m
├── cca_ssvep.m
├── generate_ssvep_references.m
├── run_offline_cca.m
├── run_online_cca.m
└── README.md
```

---

# Processing Pipeline

1. Raw EEG acquisition
2. Epoch extraction from stimulus markers
3. Bandpass filtering (4–45 Hz)
4. Reference signal generation (frequency + harmonics)
5. CCA-based frequency classification

---

# Offline vs Online Processing

## Offline Analysis

Offline processing assumes the full EEG recording is available. Each trial is analyzed independently using a fixed time window (5 seconds in this implementation).

Typical workflow:

Trial → extract EEG window → run CCA → predict target

Run offline analysis:

```matlab
results = run_offline_cca('processed/Data_Offline5.mat', ...
                          'processed/Label_Offline5.mat');
```

---

## Online Simulation

Online BCIs receive EEG samples sequentially from the acquisition system. A decision can only be made once enough samples are collected.

This repository simulates that behavior using a buffer-based approach.

Workflow:

EEG stream → buffer → CCA → predicted target → reset buffer

Current implementation:

- buffer length: **1 second**
- non-overlapping windows

Run simulation:

```matlab
predictedLabels = run_online_cca('Data.mat', 'Label.mat');
```

---

# Target Frequencies

Five SSVEP targets are used:

| Target | Frequency |
|------|------|
| 1 | 6.67 Hz |
| 2 | 7.50 Hz |
| 3 | 8.57 Hz |
| 4 | 10.00 Hz |
| 5 | 12.00 Hz |

Phase offsets:

```
[0, 0.5π, 0.5π, 0, π]
```

Reference templates include the fundamental frequency and **three harmonics**.

---

# Data Assumptions

The raw `.mat` file must contain:

```
data
```

Structure assumptions:

- EEG and marker channels stored in a matrix
- marker channel index = **26**

Trigger codes:

| Trigger | Meaning |
|------|------|
| 4 | Flicker-active period |
| 3 | Target transition |

The target label is encoded in the sample **preceding the target-transition marker**.

---

# Expected Dataset Format

After preprocessing:

```
data   → [nTrials × nChannels × nSamples]
labels → [nTrials × 1]
```

---

# How to Run

## 1. Preprocess raw data

```matlab
prepare_offline_data('path/to/raw_recording.mat', 'processed');
```

## 2. Offline classification

```matlab
results = run_offline_cca('processed/Data_Offline5.mat', ...
                          'processed/Label_Offline5.mat');
```

## 3. Online simulation

```matlab
predictedLabels = run_online_cca('Data.mat', 'Label.mat');
```
