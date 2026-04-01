# 13C_BrainMets_NABP
#### Scripts to generate the results and figures reported in the manuscript entitled "Elevated hyperpolarized [1-13C]-lactate-to-pyruvate ratio in the brain during cancer metastasis and treatment" 

##### Main_Script.R
The Main_Script.R file should be run with R. It is used to run the mixed effects linear regression model selection and final models reported in the manuscript. This R script is also used to generate the figures in Figure 2, 3 and 5 (files of which are saved in their respective directories in this repository).

##### Input Data
The data used as inputs to Main_Script.R are provided in .csv files and are described below. 

- BrainMetabolismMasterDataSheet.csv includes all average 13C signals, 13C metabolite ratios, grey matter density and volume, white matter density and volume, and brain region volume for each SLANT brain reigon (refered to as segmentation ID or SegID for short). Please note that age and sex data are removed from this dataset and therefore, the scripts for this work cannot be run unless written agreements between the corresponding author and requesting author Research Ethics Boards / Institutional Review Boards are in place.

- GTVInvolvedBrainRegions_Patient.csv reports the brain regions (SegID) that interected with gross tumor volume per patient participant.

- SegID_Volumes.csv reports the volume of each SLANT brain region.

- SegID_WMGMid.csv reports two descriptors for each brain region. First is whether the brain region is grey matter, white matter, both grey and white matter or ventricle. The second defines which lobe the brain region is apart of (frontal, parietal, occipital, temporal, insular, subcortex or white matter).


##### Notation:
- SegID = segmentation ID = SLANT brain region ID number
- PVE1 = partial volume estimate 1 = grey matter density
- PVE1_Vol = grey matter volume
- PVE2 = partial volume estimate 2 = white matter density
- PVE2_Vol = white matter volume
