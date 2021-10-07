These scripts are quite messy. If I don't absolutely have to, I won't revisit these scripts.  
- random subsets differ between surveys due to the quality of the gps-data. Therefore, I couldn't create simple scripts that woprk for each survey.  
- not everything is reproducible, in particular there seems to have been an issue with set.seed, so often the random subset cannot be reproduced  
- PS81 & PS96 have multiple files due to issues with the first subset. This is also messy.
- I tried to do too many things within a single script, making them difficult to run again.

Each file contains code that aims to:
- assign coordinates (lon/lat) to each image
- identify which images are bad/blurry
- subset dataset to a standardised number of images per transect length or per station
- crop edges of images if needed
- save images into a new directory for annotation
- save metadata in a consistent format
