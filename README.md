# ARC_Data

Main Script is:  
__ReadIn_Circumpolar_Environmental_Data.Rmd__  

All scripts in the main folder are in Rmarkdown, so each file has a _.Rmd_ and a _.html_.

All raw and derived environmental files are in two folders on Owncloud within the folder _EnvironmentalData_: _raw_ and _derived_.  
https://owncloud.imas-data-service.cloud.edu.au/index.php/s/ORxSWb6xbJRWfNI  
Currently, the easiest way is to download the folders to your local machine, keep the folder structure the same and simply edit the file-path at the top of each script that points to the "EnvironmentalData"-folder (that is if you want to generate any new files).  

## Naming convention for the files:  
*  "ReadIn_": scripts/files that read in raw or derived  datas and chnage the files into a format ready for analysis  
*  "EnvData_", "DownwardImages_", "Diatom": also self explanatory  
*  "500m","5km", etc: resolution of the raster. If none, then it's the original resolution from the raw data.  
*  "shelf": in this case the data in the raster is restricted to the depth area between 0 - 3000m  
*  "scaled": data are scaled by their mean  

NPP: net primary productivity  
waom2k: ROMS ocean model at 2km original resolution  


