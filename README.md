# EMOD_mir-184D

Single node simulation use scripts Prashanth Selvaraj, Dec 2024

Requirements: emodpy, emodpy-malaria, idmtools packages. These are available upon request from idm@gatesfoundation.org

COMPS system for HPC job management Email idm@gatesfoundation.org

Input files and executable are in the download directory. 

Files and directories:

analyzers/SummaryReportAnalyzer.py - analyze malaria epidemic data output files

analyzers/VectorGeneticsAnalyzer.py - analyze vector genetics data output files

analyzers/VectorStatsAnalyzerYearly.py - analyze vector stats data output files

download/Eradication - executable to run EMOD-Malaria 2.23

download/scehma.json - describes all the parameters available in EMOD-Malaria 2.23

input_files/single_node_demographics.json - demographics file for single node simulations to set up human population

dtk_sif.id - contains ID for docker identifier on the IDM servers

helpers.py - helper functions to set up single node simulations

idmtools.ini - contains scheduling information for the IDM servers

manifest.py - describes paths to eradication and schema file for each simulation

params.py - describes number of seeds and experiment names to use for simulations

run_mir_drive.py - creates and runs non-spatial simulations for different mir drive parameters

run_ssmt_analysis.py - main file for running analyzers