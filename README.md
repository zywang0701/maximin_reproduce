## Code

This folder contains all necessary files for real-data analysis and
simulation settings in the paper.

### Real data

-   ./code/realdata.R guides users through all the necessary steps to
    analyze the real data. Upon completion, users can replicate all
    figures and tables found in section 8 of the paper. The script takes
    more than 8 hours to finish.

### Simulations

The folders “./code/simulations” corresponds to the simulations in the
paper. These two folders contains all necessary files to generate and
analyze the simulated settings. The contents are as follows:

-   ./code/simulations/I0to6.R generates setting (I0), (I1),…, or (I6)
    and apply the proposed algorithm on it
-   ./code/simulations/I7&8.R generates setting (I7) or (I8) and apply
    the proposed algorithm on it
-   ./code/simulations/I9&10.R generates setting (I9) or (I10) and apply
    the proposed algorithhm on it
-   ./code/simulations/Setting1.R generates setting 1 and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting2.R generates setting 2 and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting3a.R generates setting 3a and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting3b.R generates setting 3b and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting4.R generates setting 4 and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting5.R generates setting 5 and apply the
    proposed algorithhm on it
-   ./code/simulations/Setting6.R generates setting 6 and apply the
    proposed algorithhm on it

Each script (setting) listed above requires more than 8 hours to finish.

### Source

The folder “./source” contains all source code. The contents are:

-   ./source/MaximinInference.R: main functions related to Maximin
    Inference
-   ./source/RDsource.R: functions to analyze real data
-   ./source/gendata.R: util functions used to generate simulated data
-   ./source/LFsource.R: bias-correction approach to handle
    high-dimensional regression
-   ./source/utils.R: other helper functions

These functions uses the following libraries:

-   glmnet version 4.1.4
-   CVXR version 1.0.11
-   intervals version 0.15.2
-   MASS version 7.3.57
