# Social_networks_kin_recognition

The “Social_networks_kin_recognition” folder contains supplemental code for Scott, Grafen & West (in review). All code is written for MATLAB R2019a.

There are three subfolders. The subfolder "Scripts_that_generate_data" contains three scripts. These scripts generate data by iterating our population genetic models (recursions). One script generates data for the “simple model” described in Supp. Info. 3. One script generates data for the weak selection version of our “island model”, which is described in the main text & Supp. Info. 4. One script generates data for the stronger selection (agent-based simulation) version of our “island model”, which is described in Supp. Info. 4d.

The subfolder "Saved_data" stores data that was generated using the scripts from the "Scripts_that_generate_data" subfolder. These data files underpin the figures in the main text and Supp. Info.

The subfolder "Scripts_that_use_data_to_generate_figures" contains multiple scripts. The scripts contain code used to generate figures featured in the main text and Supp. Info. Some of these scripts load data from the "Saved_data" folder, and others do not load data (self-contained). Code is given for each figure in the main text & Supp Info, aside from figures that: (i) are not based on data; (ii) plot simple functions given in the text.
