import xarray as xr # import the xarray package

# load waom4 depth field from sigma-layers parameters:
ds = xr.open_mfdataset(paths='E:/science/data_environmental/Circumpolar_ROMS/4km_outputs/output_yr10/ocean_avg_0001.nc')

# ds is the xarray dataset from the model output (ocean_avg_0010.nc)
ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'h', 'Vtransform']) # this line set the coordinates (not sure this is necessary in other programming languages)

# this is the calculation you’re interested:
Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho + ds.zice

# writing in other way without xarray:
#Zo_rho = (hc*s_rho + Cs_r*h) / (hc + h) # here, hc, Cs_r, hc are parameters (scalars) values for the terrain-following vertical coordinates
# s_rho is the vertical levels dimension (varying from 0 to 1)
# h is the model bathymetry (2d)
#Z_rho = zeta + (zeta + h) * Zo_rho * zice
# here zeta is the free-surface elevation (3d var)
# zice is the ice draft, needs to be accounted to get the full depth (in case the grid point sits within the ice shelf cavity).

del Zo_rho # delete tmp variable
ds.close() # close xarray dataset


