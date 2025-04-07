# Coexistence_Fluctuation
Code for first chapter of PhD thesis

Heerkou_Pre.csv: Rainfall data per month from 1901-2020

Competition.csv: Seed data

s_g_data.csv: seed soil bank survival and germination rate

Ricker model fit.R:  Output file: pars_ricker.csv; post1000.csv

Ricker model fit random effect.R: Output file: pars_ricker_re.csv

Rain_analysis.R: Visualization of rainfall pattern across years 1901-2020 【Figure S1】

Parameters visualization.R: Code for visualization of alpha and lambda 
【Figure S4】

Parameter_compare_randomeffect.R:  Code for 【Figure S2】

igr_sample_1000.R: Draw 1000 samples from the posterior distribution and obtain the possibility of the coexistence. Output file:  list.consistent_1000.rds, list.historic_1000.rds, co_compare.csv

Coexistence_consistent.R: Code for visualization of coexistence outcome of co_compare.csv with pie chart. Output file: Figure/Co_piechart/.
【Code for Figure 1】

partitioning_1000.R: Results of Ellner's partitioning for 1000 samples. Output file: par1000.csv; partitioning_1000.rds

Partitioning_visualization.R: Visualization of comparison between delta0 and delta r and Visualization of comparison between delta_alpha and delta_lambda and delta_alpha lambda
