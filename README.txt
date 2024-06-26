##Intertemporal Meditation Regulates Time Perception and Emotions—An Exploratory fNIRS Study##

Welcome to the GitHub repository for our research project. This repository contains the code and data used in our study.

##Data and Code Usage##

##Behavioral Data Analysis##

Data: Behavioral_data.xlsx
Script: data_analysis_behavior.R

- Run the script segment by segment to follow each step of our analysis.
- Generated graphs: time_pic.pdf, intertemporal_emo_pic.pdf, mindfulness_emo_pic.pdf, contrast_emo_pic.pdf
- Final figure: Fig. 3 (integrated using Adobe Illustrator)

##Hemodynamic Data Analysis##

Data: Folders session 1, session 2, and session 3
Script: data_analysis_fnirs.m

- Run in MATLAB (version R2017b) with Homer3 toolboxes (https://doi.org/10.1364/AO.48.00D280).
- Data info and exclusions: fNIRS_participant_list.xlsx
- Preprocessing summary: ProcStreamFunctionsSummary (in each session folder)
- Generated graphs: Pair2.pdf, Pair4.pdf, Pair6.pdf, Pair8.pdf, Pair9.pdf, Pair10.pdf
- Final figure: Fig. 5 (integrated using Adobe Illustrator)

##Brain Maps Visualization##

Data for plotting: fNIRS_visualization.xlsx
Script: visualization_fnirs.m

-Toolboxes: EasyTopo (please refer to: 'EasyTopo 2.0 User Manual.pdf')
-Generated data files: input_contrast1_data.mat, input_contrast12_data.mat, input_contrast123_data.mat
-Visualized maps: visualization_contrast1.eps, visualization_contrast12.eps, visualization_contrast123.eps
-Final figure: Fig. 4 (integrated using Adobe Illustrator)

##Spatial Registration Results##

Files: SpatialRegistration_BA.txt, SpatialRegistration_MNI.txt
Toolboxes: SPM8 and NIRS-SPM (MATLAB version R2013b)
Reference: https://doi.org/10.1016/j.neuroimage.2008.08.036

##Edited on June 26, 2024