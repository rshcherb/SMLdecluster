## Supervised Machine Learning (SML) declustering: Matlab computer scripts to perform declustering of earthquake catalogs

#### Directory structure
**MATLAB** folder contains the function scripts that perform the analysis.

**italy** and **southcal** folders contain the output folders where the results of the analysis is written.

#### To perform declustering:
1. run fit_etas2d_catalogs.m to fit the ETAS model to a chosen seismicity catalog.
2. run sim_etas2d_catalogs.m to simulate the ETAS model using the parameters estmated at the previous step.
3. run decluster_catalog.m to decluster the chosen catalog using one of the four methods:
4. run ML_NND_analysis.m to perform analysis using the three declustering methods
5. run nnd_analysis.m to perform NND analysis of a given catalog.

These entry scripts specify all the parameters needed to perform the corresponding tasks.

One can choose 

The scripts are written for the following publication: R. Shcherbakov, S. Kothari (2025) Earthquake declustering using supervised machine learning
