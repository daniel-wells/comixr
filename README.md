# digi-pico
Analysis of digi-pico sequencing of Ovarian Cancers

The data is from Ahmed Ahmed's lab in the WIMM, Oxford, and this analysis is supervised by Chris Yau (Department of Statistics / WTCHG, Oxford).

The aim of the project is to develop a visalisation of somatic copy number abberations from near single cell sequencing (10-20 cells) of chemotherapy resistant ovarian tumours.

### Project directory/file structure

Functions to fit the models and plot resulting components are in the files in the [R/](R/) directory.

[changepoint.R](changepoint.R) is the main analysis script, it loads the data, calculates changepoints, and then fits the model.

[GMM-EM.pdf](GMM-EM.pdf)/[tex](GMM-EM.tex) describe the mathematical model of the shared component gaussian mixtures and the expectation maximisation updates.

[mclust-mixture-model.R](mclust-mixture-model.R) fits GGMs to individual segments as well as the whole genome, which gives initialisation values for the changepoint.R analysis.

In [vignettes/](vignettes/) there are simple examples of how to use the GMM-EM functions.

[unused.segmentations.R](unused.segmentations.R) has alternative changepoint package analyses to changepoint.R