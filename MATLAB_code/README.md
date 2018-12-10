DIVERS: Decomposition of Variance Using Replicate Sampling and Spike-in Sequencing

Software requirements: MATLAB_R2016a

Installation: Only installation of MATLAB_R2016a is required. This code has been tested on MATLAB_R2016a.

Demo: The DIVERS.zip file contains the step by step workflow for performing the DIVERS variance and covariance decomposition from the original fecal OTU time series data. The DIVERS.zip file contains three folders: DIVERS_files, DIVERS_scripts and DIVERS_figures. 

1) Open MATLAB_R2016a.

2) Type into the MATLAB command line: addpath(genpath('/Path/To/.../DIVERS/‘))
Where the user must specify the full directory path containing the folder /DIVERS. 
This will add all contents of the DIVERS folder into the MATLAB path.
Please do not copy and paste. Type this command manually.


3) Open the MATLAB script, gut_table.m contained in the folder, /DIVERS/DIVERS_files . 
	3a) Specify the full directory path containing /DIVERS/DIVERS_files/ where indicated at the top of the script. Only modify ‘/Path/To/…/’ with the appropriate path.
	3b) Save the script.

4) Type into the MATLAB command line: gut_table
This runs the script gut_table.m. After several seconds, you should see an output MAT file appear in the same folder, gut_table.mat
This script simply reads in and parses the original fecal OTU table and associated metadata, and saves the output to gut_table.mat. This file will be loaded in for downstream analyses.

5) Open the MATLAB script, DIVERS_gut.m, contained in the folder, /DIVERS/DIVERS_scripts. 
	5a) Specify the full path directory containing both /DIVERS/DIVERS_files/ and /DIVERS/DIVERS_scripts/ where indicated at the top of the script. Only modify ‘/Path/To/…/’ with the appropriate path. 
	5b) Save the script.

6) Type into the MATLAB command line: DIVERS_gut
This runs the script DIVERS_gut.m, which performs the variance and covariance decomposition of fecal OTU abundances. 
Within ~60s, you should see a text file, DIVERS_gut.txt appear in folder /DIVERS/DIVERS_scripts/. This contains the output for the DIVERS variance decomposition for all individual OTUs.
You should also see an output MAT file appear in the same folder, DIVERS_gut.mat
This output can be used to re-plot all the figures in the main text (see below).

7) To plot Figures 1-3 of the main text, type into the MATLAB command line:
	7a) figure1bcd - this will plot figures 1b,c,d of the main text
	7b) figure2 - this will plot figure 2 of the main text
	7c) figure3ab - this will plot figures 3a,b of the main text
	7d) figure3cd - this plot figures 3c,d, of the main text
