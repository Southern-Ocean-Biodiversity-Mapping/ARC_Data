# ARC_Data

Main Script is:  
ReadIn_Circumpolar_Environmental_Data.Rmd

All raw and derived environmental files are in two folders on Owncloud within the folder "EnvironmentalData".  
Easiest way is to download the folders to your local machine and edit the file-path at the top of the script if you want to generate any new files.  

Naming convention for the files:
"Circumpolar_EnvData_" : all files start like this because they are environmental data.  
"500m","5km", etc: resolution of the raster. If none, then it's the original resolution from the raw data.  
"shelf": in this case the data in the raster is restricted to the depth area between 0 - 3000m  
"scaled": data are scaled by their mean  

NPP: net primary productivity  
waom2k: ROMS ocean model at 2km original resolution  


