## Simulated Data

### Settings Generation

The folder “Irregular-Settings” and “main-Settings” correspond to the
simulations in the paper. These two folders contains all necessary files
to generate the simulated settings.

The contents are as follows:

-   ./Irregular-Settings/I0to6.R: functions to generate and simulate on
    settings (I0), (I1),…, (I6)
-   ./Irregular-Settings/I7&8.R: functions to generate and simulate on
    settings (I7), (I8)
-   ./Irregular-Settings/I9&10.R: functions to generate and simulate on
    settings (I9), (I10)
-   ./Irregular-Settings/figure3.R: function to plot the figure 3 in the
    paper
-   ./main-Settings/Setting1.R
-   ./main-Settings/Setting2.R
-   ./main-Settings/Setting3a.R
-   ./main-Settings/Setting3b.R
-   ./main-Settings/Setting4.R
-   ./main-Settings/Setting5.R
-   ./main-Settings/Setting6.R

Each script (setting) listed above requires approximately 1 to 2 minutes to 
complete a single simulation. In the paper, we reported the average results 
obtained from 500 simulations. Therefore, to reproduce these results, the total 
run time could be several hours.

### Code

The folder “source” contains all source code needed for simulations. The
contents are:

-   ./MaximinInference.R: main functions related to Maximin Inference
-   ./gendata.R: util functions used to generate data
-   ./source\_functions.R: bias-correction approach to handle
    high-dimensional regression
-   ./helper\_functions.R: other functions

These functions uses the following libraries:

-   glmnet version 4.1.4
-   CVXR version 1.0.11
-   intervals version 0.15.2
-   MASS version 7.3.57

## Real Data

The folder “Real Data” encompass the steps to process real data files
and how to obtain the final results.

-   ./RealData/step-1.R: pre-process the real data files
-   ./RealData/step-2.R: Table-3 in the paper
-   ./RealData/step-3.R: Figure-6 in the paper
-   ./RealData/step-4.R: Table-4 in the paper
-   ./RealData/step-5.R: Figure-7 in the paper

The analysis of real data involves utilizing the source code located in
the folder “RealData/Maximin\_src”. While most of the code is a replica
of the source code used in simulations, we made slight adjustments to
some functions to accommodate real data.

Step-1 is time-consuming, as it requires bias-correction to be performed 
in all 46 medias, and takes more than 10 hours to complete. 
Once the results from step-1 are obtained, the remaining steps can be completed 
in just a few minutes.
