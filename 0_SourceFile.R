
if(usr == "VM"){
  usr.main.dir <- "/pvol/"
  usr.roms.dir <- "/pvol3TB/"
  usr.dropbox.dir <- "/pvol/DropboxData/"
}

if(usr =="SJ"){
  usr.main.dir <- "~/Data"
}

if(usr == "JJ"){
  usr.main.dir <- "C:/Users/jjansen/UTAS Research Dropbox/Jan Jansen/Data/"
}

## polar stereographic projection for Antarctica:
#stereo <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
stereo <- "EPSG:9354"

