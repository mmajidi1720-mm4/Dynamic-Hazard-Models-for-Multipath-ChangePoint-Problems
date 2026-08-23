# Implementation and Reproducibility Package for "Dynamic Hazard Model for Multipath Changepoint Detection"

This repository provides the complete R script necessary for **both** the practical application of the proposed Dynamic Hazard Changepoint Model (DHCM) to real-world datasets, and the reproduction of the simulation studies presented in the manuscript titled, *"Dynamic Hazard Models for Multipath ChangePoint Problems: Theory, Inference, and Application to Cognitive Decline"*.

Based on peer-review feedback, the code has been significantly updated and consolidated into a single, unified, and well-documented script (`Dynamic_Hazard_Changepoint_Model_R4.R`). 

## Repository Contents

*   **`Dynamic_Hazard_Changepoint_Model_R4.R`**: The sole R script in this repository, structured into two distinct sections:
    *   **Part I: Practical Implementation Guide:** Contains the generalized `fit_dhcm` function, which dynamically accommodates any number of covariates and estimates the *full* parameter vector (including all baseline and post-changepoint distributional parameters). It also includes `predict_changepoints` for individual-level estimation and a minimal, self-contained example demonstrating how to format data and fit the model.
    *   **Part II: Simulation Reproducibility:** Contains the routines required to reproduce the Monte Carlo simulations (Section 4 of the manuscript), comparing the proposed DHCM against the Constant Hazard Model and the PELT algorithm.

## System and Software Requirements

*   **R Version:** Developed and tested on R 4.2.0 and later.
*   **R Packages:**
    *   `stats` (included in base R)
    *   `parallel` (included in base R - used for multi-core simulation processing)
    *   `changepoint` (optional, for PELT comparison in Part II. The script includes a fallback CUSUM method if this package is not installed).

You can install the optional `changepoint` package with:
```R
install.packages("changepoint")
```

## How to Use the Code

### A. Practical Application to Your Own Data (Part I)
If you wish to apply the DHCM to your own longitudinal dataset:
1. Open `Dynamic_Hazard_Changepoint_Model_R4.R`.
2. Run the functions defined in **Part I**.
3. Follow the "Minimal Self-Contained Example" provided at the end of Part I. It demonstrates how to structure your design matrix (`X_cov`), run `fit_dhcm()`, and extract individual change-point estimates ($\hat{\tau}_i$) using `predict_changepoints()`.

### B. Reproducing Manuscript Simulations (Part II)
To reproduce the exact simulation results reported in the paper:
1.  **Download and Open**: Download `Dynamic_Hazard_Changepoint_Model_R4.R` and set your R working directory to its location.
2.  **Choose the Simulation Mode**: In **Part II** of the script, look for the `FAST_MODE` toggle:
    *   **Fast Test Run**: To quickly verify the code runs without errors, set `FAST_MODE <- TRUE`. This runs a small number of replications (`R_sim = 50`) and completes in minutes.
    *   **Full Manuscript Replication**: To reproduce the full tables from the paper, set `FAST_MODE <- FALSE`. This runs 1,000 replications using parallel processing. *Note: This will take a significant amount of time depending on your CPU.*
3.  **Execute**: Run the script using the following command in your console:
```R
source('Dynamic_Hazard_Changepoint_Model_R4.R')
```
