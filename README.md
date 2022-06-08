# Social_networks_kin_recognition

The “Social_networks_kin_recognition” folder contains supplemental code for "Multiple social encounters can eliminate Crozier’s paradox and stabilise genetic kin recognition" (Scott, Grafen & West; Nature Communications; 2022). All code is written for MATLAB R2019a. There are three subfolders.

The subfolder "Scripts_that_generate_data" contains three scripts. The "Generate_data_for_island_model_weak_selection" script generates data for the weak-selection mathematical model. The "Generate_data_for_agent_based_simulation" script generates data for the agent-based simulation model, which accounts for stronger selection and finite populations. The "Generate_data_for_balancing_selection_finite_pop" script generates data for the version of the agent-based simulation model where there is no tag mutation, in which balancing selection is examined via tag fixation times.

The subfolder "Saved_data" stores data that was generated using the scripts from the "Scripts_that_generate_data" subfolder. These data files underpin the figures in the main text and Supp. Info.

The subfolder "Scripts_that_use_data_to_generate_figures" contains multiple scripts. The scripts contain code used to generate figures featured in the main text and Supp. Info. Some of these scripts load data from the "Saved_data" folder, and others do not load data (i.e. they are self-contained). Code is given for each figure in the main text & Supp Info, aside from figures that: (i) are not based on data; (ii) plot simple functions given in the text.
