# Connectivity Model

This repository contains code for the simulation and analysis of the research project:
**"Connectivity indicators can predict population persistence in river networks"**.

## Overview

The project provides tools and simulation scripts to analyze how structural and functional
connectivity in river networks influences the long-term persistence of populations. 
It focuses on quantifying connectivity indicators and exploring their predictive power
for metapopulation persistence under varying network structures and dispersal scenarios.

## Features

- **Simulation of River Network Connectivity:**  
  Scripts for generating and analyzing both linear and binary river network topologies,
  including symmetric and asymmetric (directional) dispersal.

- **Computation of Connectivity Indicators:**  
  Calculates various indicators such as DCI (Degree of Connectivity Index), 
  reproductive value, and local steady state occupancy for network reaches.

- **Barrier Analysis:**  
  Tools to simulate the removal of barriers and assess their impact on global
  and local connectivity indicators.

## Structure

- `scripts/`: Main R scripts for simulation and analysis:
  - `simData.R`: Generates and manages simulation data.
  - `global.R`: Analyzes global indicators across network configurations.
  - `local1.R` and `local2.R`: Analyze local indicators and effects of barrier removals.
  - `fun_con.R`: Core functions for connectivity and metapopulation calculations.
  - `netwSimInputs.R`: Parameter definitions and simulation input setup.

- `data/`: Expected location for input and output data files (not included in the repository).

## Requirements

- R (recommended >= 4.0)

Additional dependencies may be installed as needed by the scripts.

## Usage

1. Clone the repository:
   ```sh
   git clone https://github.com/aligharouni/connectivity_model.git
   ```

2. Install required R packages (see script headers for details).

3. Run the scripts in `scripts/` for data generation, analysis, and visualization. 
Adjust simulation parameters in `netwSimInputs.R` as necessary.

## Citation

If you use this code or its results in your research, 
please cite the corresponding publication (project title above).

---

