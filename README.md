# Venus_O2

Code and datasets required to reproduce figures in main text of Warren & Kite 2023 "Narrow range of early habitable Venus scenarios permitted by modelling of oxygen loss and radiogenic argon degassing." https://www.pnas.org/doi/full/10.1073/pnas.2209751120

To download data used to produce publication figures, go to: https://zenodo.org/record/7416204#.Y5UVEXbMK3A

After downloading, code will need to be modified to find files in chosen directories for running the full model (also requires installation of VolcGases from github.com/Nicholaswogan/VolcGases) and files for re-creating plots. All .zip files contain data used to create figures in paper. 


Instructions to reproduce figures:

1. To reproduce full dataset, download all files excluding .zip files and install VolcGases. Modify "modular_functions_clean_redox.py" to match VolcGases installation.

  For runaway greenhouse runs: Ensure "modular_functions_clean_redox.py" line 512 is commented out. Run "adding_dissolution_clean.py" followed by "ext1line.py" to pre-run runaway greenhouse model and save output to read into full model (this speeds up running the entire suite of "melting" models, but is not strictly necessary). Next, use "gridsearch_all_redox.csv" to set model parameters, then run "input_clean_redox.py" to initiate model. 

  For runs without runaway greenhouse melting: uncomment "modular_functions_clean_redox.py" line 512. Use "gridsearch_all_redox.csv" to set model parameters, then run "input_clean_redox.py" to initiate model. 

2. To reproduce Figures 2, 3, S1, and S2, either download all zip files beginning with "Fig2_" and "Fig3_", extract all files to single location, or save all new model output into a single directory and use "plot_together.py" to generate figures (instructions contained within script). Use same script to reproduce Figures S9 and S10 with data in "t_atm_sensitivity.zip" and "bH_sensitivity.zip". 

3. To reproduce Figures 4 and 5, If using new model runs to generate plots, first run "gen_data_forplots.py". Alternatively, download:

  - e_statistics_nomelt_CO_FMQ0.npz 
  - e_wd_statistics_melt_CO_FMQ0.npz
  - 40_Ar_hab_dlith_KU_mix_52522.csv
  - 40_Ar_hab_dlith_KU_nomix_52522.csv

  Then, run maintext_plotting_new.py.

4. To reproduce Figures S5 to S8, run "Ar_plots.py". Requires: "serpent_dehyd2.csv" and "eclogite_transition.csv".
