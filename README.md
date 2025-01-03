## Supervised Machine Learning (SML) declustering: Computer scripts to perform declustering of earthquake catalogs

The scripts were written for the following publication: R. Shcherbakov, S. Kothari (2025) Earthquake declustering using supervised machine learning, submitted to BSSA.

### Directory structure
**MATLAB/** folder contains the Matlab scripts that perform specific tasks. The following subfolderes contain functions:
- **decluster/** to perform declustering.
- **etas2d/** to fit and simulate forward the 2D ETAS model.
- **fm/** to analyze the frequency-magnitude statistics.
- **ml_nnd/** to perform the SML declustering.
- **nnd/** to perfomr NND analysis.

**italy/**, **southcal/** folders contain the output subfolders where the results of the analysis is written for the specified region. It is possible to add other regions. These folders have the following subfolders:
- **catalog/** must contain an earthquake catalog for the specified region. It must be a text file with 10 columns: \[year month day hour min sec latitue longitude magnitude depth\]
- **decluster/** 
- **etas_fit/**
- **etas_test_cat/** 
- **etas_train_cat/**
- **ml_nnd/**
- **nnd/** 

**etas2d_parameters.m** script contains information about the geographical boundaries of the full and target regions, time intervals, initial conditions, and other parameters that are specific for each region.

### To perform declustering:
1. run **fit_etas2d_catalogs.m** to fit the ETAS model to a chosen seismicity catalog.
2. run **sim_etas2d_catalogs.m** to simulate the ETAS model using the parameters estmated at the previous step.
3. run **decluster_catalog.m** to decluster the chosen catalog using one of the four methods:
4. run **ML_NND_analysis.m** to perform analysis using the three declustering methods
5. run **nnd_analysis.m** to perform NND analysis of a given catalog.

These entry scripts specify all the initial parameters needed to perform the corresponding tasks. **Model** structure is used to pass the information into most functions. Model.sRegion specifies the region to analize and must match the corresponding folder name. 

