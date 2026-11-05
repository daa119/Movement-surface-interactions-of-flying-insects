# Movement-surface interactions of flying insects in three-dimensional space

## Overview

This repository contains MATLAB simulation codes accompanying the research manuscript:

**Movement-surface interactions of flying insects in three-dimensional space**

**Authors:** Danish A. Ahmed, Ross N. Cuthbert, Isadora E. Fluck, Frances M. Hawkes and Gala Garrod.

The study develops a three-dimensional stochastic modelling framework to investigate how insect flight behaviour influences encounters and capture by finite planar trapping surfaces.

Individual insect movement is simulated using correlated random walks (CRWs) and biased correlated random walks (BCRWs), incorporating directional persistence, directional attraction towards a trapping panel, geometric interception, repeated panel contacts, probabilistic retention and edge leakage.

The framework examines how movement behaviour, panel geometry and retention jointly influence trapping outcomes, with potential applications in insect surveillance, pest management and vector monitoring.

## MATLAB simulation codes

The repository contains seven MATLAB scripts used to generate Figures 1–8 of the manuscript.

| File | Description |
|---|---|
| `fig_1_2.m` | Three-dimensional trajectories, movement configurations and corresponding two-dimensional projections (Figures 1 and 2). |
| `fig_3.m` | Encounter and capture probabilities as functions of panel size relative to movement persistence length (Figure 3). |
| `fig_4.m` | Contact frequency and repeated panel interactions under different bias strengths and retention probabilities (Figure 4). |
| `fig_5.m` | Capture probability heat maps examining the combined effects of directional persistence, directional bias and retention (Figure 5). |
| `fig_6.m` | Capture probability as a function of retention probability under different movement regimes (Figure 6). |
| `fig_7.m` | Mean time to first panel encounter as a function of initial distance, directional bias and persistence (Figure 7). |
| `fig_8.m` | Edge leakage probability as a function of directional persistence, directional bias and panel size (Figure 8). |

## Model parameters

The simulations use the following baseline movement parameters:

- Mean flight speed: 25 cm/s.
- Time step: 0.02 s.
- Mean step length: 0.5 cm.
- Step-length distribution: Maxwell distribution.
- Directional updating: von Mises–Fisher distribution.
- Baseline square panel side length: 10 cm.
- Number of trajectories: 1,000 per parameter combination for Figures 3–8.
- Maximum simulation steps: 2,500 for Figures 3–8.

Directional persistence, directional bias, panel dimensions and retention probabilities are varied according to the analyses presented in the manuscript.

## Software and reproducibility

The simulations were developed using **MATLAB R2016a**.

Each MATLAB program contains the relevant model functions, parameter definitions and plotting instructions.

To reproduce the figures:

1. Download the MATLAB files from this repository.
2. Open MATLAB and navigate to the directory containing the files.
3. Run the corresponding figure-generating function.
4. The simulations will generate the figures using the parameters specified in each file.

The scripts initialise random number generators with specified seeds to support reproducibility.

Simulation settings can be modified within each script to explore alternative movement behaviours, panel sizes and trapping conditions.

## Citation

If you use these codes in your research, please cite the accompanying manuscript:

Ahmed, D. A., Cuthbert, R. N., Fluck, I. E., Hawkes, F. M., and Garrod, G. *Movement-surface interactions of flying insects in three-dimensional space*.

A permanent DOI for the archived code will be added following deposition in Zenodo.

## Contact

**Danish A. Ahmed**

Department of Mathematics and Natural Sciences

Gulf University for Science and Technology, Kuwait

Email: Ahmed.D@gust.edu.kw
