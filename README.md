supp data for ...

# data
- `responses&predictions.xls` contains behavioral responses and cross-validated model predictions for each syllable pair and stimulus sequence tested (see stimulusSequences.xls)
- `stimulusSequences.xls` - the sheet "stimuluspattern" contains all 32 unique sequences of subunits tested as sequences of 0s and 1s. The sheet "subunitmapping" lists which subunit is encoded as 0 and as 1 for each pair of subunits tested.
- `modelParameters.xls` contains the mean/std/sem over the 218 cross-validaton runs for each model parameter

# code
The demo script `fit_driftdiffusionmodel.m` loads data and fits a drift-diffusion model using a genetic algorithm.
