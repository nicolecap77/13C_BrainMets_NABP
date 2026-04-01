# 13C_BrainMets_NABP
#### Scripts to generate the results and figures reported in the manuscript entitled "Elevated hyperpolarized [1-13C]-lactate-to-pyruvate ratio in the brain during cancer metastasis and treatment" 

#### Main_Script.R
The Main_Script.R file should be run with R. It is used to run the mixed effects linear regression model selection and final models reported in the manuscript. This R script is also used to generate the figures in Figure 2, 3 and 5 (files of which are saved in their respective directories in this repository).

#### Input Data
The data used as inputs to Main_Script.R are provided in .csv files and are described below. Please put all of these csv files and the custom functions (See Custom Functions section next) in the same directory as Main_Script.R

- BrainMetabolismMasterDataSheet.csv includes all average 13C signals, 13C metabolite ratios, grey matter density and volume, white matter density and volume, and brain region volume for each SLANT brain reigon (refered to as segmentation ID or SegID for short). Please note that age and sex data are removed from this dataset to protect the personal health information of participants and therefore, the scripts for this work cannot be run unless written agreements between the corresponding author and requesting author Research Ethics Boards / Institutional Review Boards are in place. Please contact the corresponding author Dr. Charles H Cunningham (charles.cunningham@utoronto.ca) to request a data agreement.

- GTVInvolvedBrainRegions_Patient.csv reports the brain regions (SegID) that interected with gross tumor volume per patient participant.

- SegID_Volumes.csv reports the volume of each SLANT brain region.

- SegID_WMGMid.csv reports two descriptors for each brain region. First is whether the brain region is grey matter, white matter, both grey and white matter or ventricle. The second defines which lobe the brain region is apart of (frontal, parietal, occipital, temporal, insular, subcortex or white matter).

#### Custom Functions
Two custom functions were made and are required to run Main_Script.R

- generate_regression_formulas.R makes all possible regression formula combinations up to and including two-way interaction terms, given the dependent and independent variables as input.

- pval_to_asterisk.R assigns different p-value levels to different number of asterisks for plots.

#### Regression Results
A directory called Regression_Results stores the results of model selection, final model results, estimates marginal means and post hoc tests (where applicable) for each dependent variable explored in the manuscript. There are 4-5 .csv files in each of the directories according.

- DependentVariable_ModelSelection_AIC.csv reports each model tested in the mixed effects linear regression model selection process, with AIC reported for each model. Models are sorted in ascending AIC order.
  
- DependentVariable_MELR_Summary.csv reports the mixed effects linear regression (MELR) results of the optimal model selected by in the model selection process.

- DependentVariable_MELR_EMM.csv reports the estimated marginal mean value of the dependent variable between patients and controls after accounting for other independent variables using the mixed effects linear regression chosen in DependentVariable_MELR_Summary.csv. 'emmean' columns reflect dependent variables that have been log transformed while 'response' columns reflect dependent variables that have not been transformed.

-  DependentVariable_MELR_EMM_t-test.csv reports the estimated difference in dependent variable between patients and controls in the emmean or response column. 'ratio' columns report ratios between patients and controls (for log transformed dependent variables) while 'estimate' columns report differences (for non-transformed dependent variables).

-  DependentVariable_lobes_tukey.csv reports tukey HSD tests between estimated marginal mean of dependent variable grouped by brain lobe.

#### Figures:
Figures are provided in .tiff and .eps format

#### A Note on Notation:
- SegID = segmentation ID = SLANT brain region ID number
- PVE0 = partial volume estimate 0 = cerebral spinal fluid PVE
- PVE0_Vol = volume of cerebral spinal fluid
- PVE1 = partial volume estimate 1 = grey matter PVE
- PVE1_Vol = grey matter volume
- PVE2 = partial volume estimate 2 = white matter PVE
- PVE2_Vol = white matter volume
