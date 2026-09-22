Data and code to replicate the statistical analyses and figures of the article “Gaze and perception covary with reading direction”, https://doi.org/10.64898/2026.03.04.707660
All analyses described in the manuscript are included in the **Analysis** folder below, in the form of four scripts plotting figures 2 to 5 in the manuscript. These scripts run a full analysis of the single-trial data provided in the **Data** folder. 

### Analysis folder
contains R scripts used to run analyses and produce corresponding manuscript figures
 - figures[2-5]_*_git.R: compute the statistical models for and plot those respective figures, calling on the single-trial files in the **Data** folder
 - **eyetracking_preprocessing** folder contains two R scripts used to pre-process the raw eyetracking data, resulting in eyedat_clean.csv provided in the **Data** folder. Raw eye-tracking data are not provided but can be requested from the authors.
 - **figures** folder contains the figures produced by the above R scripts
 - **PNG** folder contains image files used as background in figures 3, 4
 - ovalCoordinates.csv: used for some plots

### Data folder
contains compiled single-trial data files
- biasdat_clean.csv: single-trial behaviour data for all subjects
- langdata_clean.csv: language and handedness data for all subjects
- eyedat_clean.csv: single-trial eyetracking + behaviour data for 115 eye-tracking subjects

### Materials folder
contains the study materials
 - Chimeric Face_Eye Tracker.py is for the experimental task with eyetracking
 - **online_Pavlovia** folder contains the code for the online task (no eyetracking)
 - **Final Face Stimuli** folder includes a subset of the face stimuli called by the experimental scripts. The full set of faces is not covered by the open license, but may be requested from the authors.
