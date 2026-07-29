# Reproducibility Package for "A Dynamic Hazard Model for Multipath Changepoint Detection"

This repository provides the R script necessary to reproduce the simulation studies presented in **Section 5** of the manuscript titled, *"Dynamic Hazard Models for Multipath ChangePoint Problems: Theory, Inference, and Application to Cognitive Decline"*.

The code is provided as a single, self-contained script (`Dynamic_Hazard_Changepoint_Model.R`) to ensure a straightforward and transparent verification process.

## Repository Contents

- **`Dynamic_Hazard_Changepoint_Model.R`**: The sole R script in this repository. It contains all functions and procedures required to generate simulated data, fit the proposed Dynamic Hazard Model, fit the comparative Constant Hazard Model, run the PELT algorithm, and generate the summary tables (Tables 1, 2, and 3) as seen in the manuscript.

## System and Software Requirements

- **R Version**: Developed and tested on R `4.2.0` and later.
- **R Packages**:
    - `stats` (included in base R)
    - `parallel` (included in base R)
    - `changepoint` (optional, for PELT comparison. The script includes a fallback CUSUM method if this package is not installed).

You can install the optional `changepoint` package with:
```R
install.packages("changepoint")

## How to Reproduce the Results

1.  **Download the Script**: Download the `Dynamic_Hazard_Changepoint_Model.R` file from this repository to your local machine.

2.  **Set the Working Directory**: Open R or RStudio and set your working directory to the location where you saved the script.

3.  **Choose the Simulation Mode**:
- **Fast Test Run**: To quickly verify that the code runs without errors, set `FAST_MODE <- TRUE` at the top of the script. This will run a small number of replications (`R_sim = 50`) and should complete in a few minutes.
- **Full Manuscript Replication**: To reproduce the exact results reported in the paper, set `FAST_MODE <- FALSE`. This will run the full number of replications (`R_sim = 1000`) using parallel processing on multiple CPU cores. This may take a significant amount of time depending on your hardware.

4.  **Execute the Script**: Run the entire script in your R console:

```R
    source('simulation_main.R')
