# README

This repository contains data and code for analyzing population dynamics and species coexistence under varying rainfall scenarios. Below is an overview of the files and their purposes.

## Data Files

1. **Heerkou_Pre.csv**  
   - Contains monthly rainfall data spanning from 1901 to 2020.

2. **Competition.csv**  
   - Contains seed count data derived from competition experiments.

3. **s_g_data.csv**  
   - Includes seed soil bank survival and germination rates.

## Script Files

### 1. **Ricker model fit.R**  
   - Fits population model parameters using the Ricker model.  
   - **Output:**  
     - `post_1000.rds`  
     - `pars_ricker.csv`  

### 2. **Ricker model fit random effect.R**  
   - Fits population model parameters with "pot" as a random effect.  
   - **Output:**  
     - `pars_ricker_re.csv`  

### 3. **Parameters visualization.R**  
   - Visualizes the estimated parameters.  
   - **Output:**  
     - `Figure S4`  

### 4. **Parameter_compare_randomeffect.R**  
   - Compares population model parameters with and without "pot" as a random effect.  
   - **Output:**  
     - `Figure S2`  

### 5. **igr_constant_1000.R**  
   - Simulates population dynamics and calculates the invasion growth rate (IGR) under four constant rainfall scenarios.  
   - **Output:**  
     - `list.consistent_1000.rds`  

### 6. **igr_partitioning_1000.R**  
   - Partitions the invasion growth rate into different components.  
   - **Output:**  
     - `partitioning_1000.rds`  
     - `par1000.csv`  

### 7. **Visualization_6 scenarios.R**  
   - Visualizes species coexistence outcomes under six different scenarios:5 constant rainfall scenarios (consistent wet, post-flood dry, pre-flood dry, consistent dry, average) and a variable interannual rainfall scenario .  
   - **Output:**  
     - `Figure 1`  
     - `Figure S6`  

### 8. **Partitioning_visualization.R**  
   - Visualizes the results of the invasion growth rate partitioning.  
   - **Output:**  
     - `Figure 2`  
     - `Figure S5`  
     - `Figure S7`  

### 9. **Possibility change.R**  
   - Visualizes how different coexistence mechanisms contribute to the coexistence probability.  
   - **Output:**  
     - `Figure 3`  

