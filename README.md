# ARC_Data

This respository contains scripts that setup the biological and environmental data for analysis. 

/ASAID_Classification...: Contains the scripts to create the classification catalogs published under https://doi.org/10.25959/J4EE-9A56
/FAM:  Contains the scripts for running and analysing food-availability models
/annotation:  scripts to prepare image annotations from the files available under https://doi.org/10.25959/2R0E-PB18
/prep_environment:  loading all environmental data from its source, matching extents and projection and resolution, scaling
/prep_image:  legacy scripts originally used to prepare image metadata before publication on Squidle+
/legacy_scripts:  scipts developed but not updated to suit new pipeline
/old:  scratch folder for old scripts to keep just in case

Raw environmental files are available from their source online.
Derived environmental files are currently on the teams private dropbox, and the relevant files will be published with the analysis paper.

## Naming convention for the files:  
*  "1_...":  Sequence of scripts to follow the pipeline
*  "EnvPrep_": scripts/files that read in raw or derived environmental data and chnage the files into a format ready for analysis  
*  "AnnData_": scripts for reading and preparing image annotation files
