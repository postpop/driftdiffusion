Data and code used for estimating drift-diffusion models to behavioral responses (doi://).

# data
- `data/responses&predictions.xls` contains behavioral responses and cross-validated model predictions for each syllable pair and stimulus sequence tested (see stimulusSequences.xls).
- `data/stimulusSequences.xls` - the sheet "stimuluspattern" contains all 32 unique sequences of subunits tested as sequences of 0s and 1s. The sheet "subunitmapping" lists which subunit is encoded as 0 and as 1 for each pair of subunits tested.
- `data/modelParameters.xls` contains the mean/std/sem over the 218 cross-validaton runs for each model parameter.

# code
The demo script `code/fit_driftdiffusionmodel.m` loads data and fits a drift-diffusion model using a genetic algorithm.

# publications
data and code used in:
>Jan Clemens, Jennifer Aufderheide (eq. contribution), Bernhard Ronacher, Evaluation of sensory cues for mating decisions in grasshoppers reveals influence of sexual selection. 2017, submitted

parts of the data (behavioral response to gap stimuli) and a previous of the code were used in:
> Jan Clemens, Stefanie Krämer (eq. contribution), Bernhard Ronacher, Asymmetrical integration of evidence during courtship decisions in grasshoppers. 2014, PNAS, 111(46):16562–16567

