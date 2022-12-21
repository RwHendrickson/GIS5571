### This script runs sensitivity analysis on the traffic dispersion parameters.
### Specifically, it will create different AADT_exposure rasters (10 meter resolution)
### by varying the sigma within the Gaussian kernel used in convolution (constant 2km radius - want more but not enough memory)

# Import libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

# System

import os
import time

# Analysis

import numpy as np # Numerical Analysis
import pandas as pd # Data Mgmt
import geopandas as gpd # Spatial Data Mgmt
import rasterio # Rasters
from rasterio.transform import Affine # Affine transformations
# Dispersion
from scipy.ndimage import convolve # Convolution <- this is cool
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.convolve.html
# Zonal Stats
from rasterstats import zonal_stats
# Silence some warnings
import warnings


# Definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

starttime = time.time()

def Save_array_to_geotiff_template(array, template_path, save_path):
    '''Saves a numpy array into a geotiff with the same CRS as the template.
    '''            

    # Get metadata from template
    rst = rasterio.open(template_path) # Open template
    meta = rst.meta.copy() # Copy template metadata
    # meta.update(compress='lzw') # Good for integers/categorical rasters
    rst.close()

    with rasterio.open(save_path, 'w+', **meta) as out: # Burn features into raster
        out.write_band(1, array)
        
# Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

# Load raster template for plotting

raster_path = os.path.join(os.getcwd(), '2_Model_Pollutant_Exposure', 'template.npy')

raster = np.load(raster_path)

# Get metadata 
template_path = os.path.join(os.getcwd(), '2_Model_Pollutant_Exposure', 'template.tif')
with rasterio.open(template_path) as rast:
    meta = rast.meta.copy() 

# Minneapolis Boundary

mpls_path = os.path.join(os.getcwd(), '1_Data_IO', 'Data', 'mpls_boundary.geojson')
mpls = gpd.read_file(mpls_path)

# rasterized aadt

aadt_path = os.path.join(os.getcwd(), '2_Model_Pollutant_Exposure', 'Rasterized_Sources', 'rasterized_aadt.tif')

aadt_rast = rasterio.open(aadt_path) # Open

aadt_band = aadt_rast.read(1)

aadt_rast.close()


# Sensitivity of Sigma ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dispersion is modeled through convolution with a Gaussian Kernel
# Wondering about Convolution?
# https://en.wikipedia.org/wiki/Kernel_(image_processing)#Convolution
# Geographic Context
# https://gis.stackexchange.com/questions/373648/seeking-arcgis-focal-statistics-alternative-in-open-source-python


# Convolution 

# Performing dispersion with an un-normalized gaussian kernel function

# Define weights kernel

def gkern(l=5, sig=1.):
    """\
    Creates gaussian kernel with side length `l` and a sigma of `sig`
    """
    ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(sig))
    kernel = np.outer(gauss, gauss)
    return kernel #/ np.sum(kernel) <- this is for standardizing.. like a focal sum

### Iterate through sigmas and convolve

print('Beginning Dispersion (Convolution) iterations')

# Initialize
sigmas = np.linspace(1, 15, 29)
print('Sigmas used: \n', sigmas)

for sigma in sigmas:
    
    # Get kernel
    kernel = gkern(200, sigma) # Gaussian kernel with 2 km diameter and sigma between 0.5 and 15
    
    # Disperse
    
    dispersion_rast = convolve(aadt_band, kernel)
    
#     # Standardize (Not doing this)
    
#     # Get indices greater with values greater/equal & less than 1 (for log)
#     small_indices = dispersion_rast<1
#     big_indices = dispersion_rast>=1
#     # Do the transform
#     std_rast = dispersion_rast.copy()
#     std_rast[small_indices] = 0
#     std_rast[big_indices] = np.log10(std_rast[big_indices])
#     # Normalize
#     std_rast[big_indices] = std_rast[big_indices]/std_rast[big_indices].max()
    
    # Normalize
      
    dispersion_rast = dispersion_rast/dispersion_rast.max()
    
    # Save
    
    savename = str(int(sigma*10)/10) + 'sig_AADT_exposure.tif'
    savepath = os.path.join(os.getcwd(),'7_Sensitivity_Analysis', 'Modeled_Dispersion', savename)
    Save_array_to_geotiff_template(dispersion_rast, template_path, savepath)
    
print('Dispersion Finished')

# Create Hazard Index Surface ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Combine with current industrial particulate matter 2.5 
### Using weights between 40 and 60 percent

print('Beginning hazard index surface creation')

# Initialize

# Get Current Industrial pm2.5 emissions (I_2020) dispersion surface
path_to_layers = os.path.join(os.getcwd(),'2_Model_Pollutant_Exposure', 'Modeled_Exposure')
path_to_I_2020 = os.path.join(path_to_layers, 'PM25Primary_current_exposure.tif')
with rasterio.open(path_to_I_2020) as rast:
    I_2020 = rast.read(1)
    
# Get Current AADT (T_2022) Dispersion Surface Filenames
# "Different sigma in gaussian kernel of 5km radius"

path_to_aadts = os.path.join(os.getcwd(),'7_Sensitivity_Analysis', 'Modeled_Dispersion')
aadt_filenames = os.listdir(path_to_aadts)
    
# Iterate through AADT Dispersion Surfaces (T_2022)
print('Weights used: \n', np.linspace(0.4, 0.6, num=5))

for aadt_filename in aadt_filenames:
    
    # Load 
    path_to_T_2022 = os.path.join(path_to_aadts, aadt_filename)
    with rasterio.open(path_to_T_2022) as rast:
        T_2022 = rast.read(1)
    
    # Iterate through weights
    
    for I_percent in np.linspace(0.4, 0.6, num=5): # Between 0.4 & 0.6 by 0.05

        # Get weights
        weights_name = str(int(I_percent*100)) + 'I-' + str(int((1-I_percent)*100)) + 'T'
        weights_test = np.array([I_percent,1-I_percent])
        
        # Calculate linear combination
        hazard_index_surface = I_percent * I_2020 + (1-I_percent) * T_2022
        
        # Save
        sigma_string = aadt_filename.split('_')[0]
        savename = sigma_string + '_' + weights_name + '_HazardIndex.tif'
        savepath = os.path.join(os.getcwd(),'7_Sensitivity_Analysis', 'Hazard_Indices', savename)
        Save_array_to_geotiff_template(hazard_index_surface, template_path, savepath)

print('Finished with Hazard Index Surface Creation')

# Compare Zonal Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# They've been saved in a format that is incompatible with rasterstats zonal stats...
# Here's the correction:

# Flip the raster
# From: https://github.com/perrygeo/python-rasterstats/issues/98

def flipud(raster, affine):
    raster = np.flipud(raster)
    affine = Affine(
        affine.a,
        affine.b,
        affine.c,
        affine.d,
        -1 * affine.e,
        affine.f + (affine.e * (raster.shape[0] - 1)),
    )
    return raster, affine

# Initialize

print('Beginning Hazard Index Zonal Statistics')

# Silence warnings for zonal_stats
warnings.filterwarnings('ignore')
# Zones Path (HOLC)
path_to_zones = os.path.join(os.getcwd(), '1_Data_IO', 'Data', 'holc.geojson')
# Load Hazard Index Surfaces
indices_path = os.path.join(os.getcwd(), '7_Sensitivity_Analysis', 'Hazard_Indices')
index_filenames = os.listdir(indices_path)

zonal_stat_gdfs = {} # Storage for all zonal stats geodataframes

# Iterate through Surfaces
for index_filename in index_filenames:
    
    # Load hazard index surface
    index_filepath = os.path.join(indices_path, index_filename)
    with rasterio.open(index_filepath) as rast:
        Hazard_Index = rast.read(1)

    # Flip the index...    
    affine = meta['transform']
    Hazard_Index_flipped, new_affine = flipud(Hazard_Index, affine) # flip the raster
    
    # Do Zonal Stats
    zs = zonal_stats(path_to_zones, Hazard_Index_flipped,
                     affine = new_affine, geojson_out=True) # Get zonal stats
    
    # Save into gdf of stats
    zonal_stats_gdf = gpd.GeoDataFrame.from_features(zs, crs = 'EPSG:26915') # Save as geodataframe
    savename = index_filename[:-4] + '_Zonal_Stats.geojson'
    savepath = os.path.join(os.getcwd(), '7_Sensitivity_Analysis', 'Zonal_Stats', savename)
    
    zonal_stats_gdf.to_file(savepath)
    
print('Done with Zonal Statisitics')
print('Total Time: ', starttime - time.time())