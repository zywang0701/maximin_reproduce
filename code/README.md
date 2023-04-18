## Code

This folder contains all necessary files for real-data analysis and
simulation settings in the paper.

### Real data

-   ./code/realdata/realdata.R guides users through all the necessary steps to
    analyze the real data. Upon completion, users can replicate all
    figures and tables found in section 8 of the paper. 

The initial section of data pre-processing in the script necessitates more than 
8 hours to complete, while the subsequent sections are completed quickly.

### Simulations

The folders “./code/simulations” corresponds to the simulations in the
paper. These two folders contains all necessary files to generate and
analyze the simulated settings. The contents are as follows:

-   ./code/simulations/I0to6.R generates setting (I0), (I1),…, or (I6)
    and applies the proposed algorithm to the generated data
-   ./code/simulations/I7&8.R generates setting (I7) or (I8) and applies
    the proposed algorithm to the generated data
-   ./code/simulations/I9&10.R generates setting (I9) or (I10) and applies
    the proposed algorithm to the generated data
-   ./code/simulations/Setting1.R generates setting 1 and applies the
    proposed algorithm to the generated data
-   ./code/simulations/Setting2.R generates setting 2 and applies the
    proposed algorithm to the generated data
-   ./code/simulations/Setting3a.R generates setting 3a and apply the
    proposed algorithm to the generated data
-   ./code/simulations/Setting3b.R generates setting 3b and applies the
    proposed algorithm to the generated data
-   ./code/simulations/Setting4.R generates setting 4 and applies the
    proposed algorithm to the generated data
-   ./code/simulations/Setting5.R generates setting 5 and applies the
    proposed algorithm to the generated data
-   ./code/simulations/Setting6.R generates setting 6 and applies the
    proposed algorithm to the generated data
-   ./code/simulations/figs&tabs.R contains functions to reproduce each
	table and figure

With the exception of the final script, each listed script (setting) necessitates 
over 8 hours to complete, as it entails a total of 500 simulations. Users 
may choose to run 500 simulations concurrently. Once the corresponding 
RData file is generated, the last script is utilized to replicate the 
figure or table, which can be completed in a matter of seconds.

### Source

The folder “./source” contains all source code. The contents are:

-   ./code/source/MaximinInference.R: main functions related to Maximin
    Inference
-   ./code/source/RDsource.R: functions to analyze real data
-   ./code/source/gendata.R: functions used to generate simulated data
-   ./code/source/LFsource.R: bias-correction methods for handling
    high-dimensional regression
-   ./code/source/utils.R: other helper functions

These functions use the following libraries:

-   glmnet version 4.1.4
-   CVXR version 1.0.11
-   intervals version 0.15.2
-   MASS version 7.3.57
