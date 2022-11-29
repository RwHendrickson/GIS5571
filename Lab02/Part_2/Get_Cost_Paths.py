### Didn't end up taking this route in the end. 
### This is a FOSS way
### Used arcpy instead - see 6_Costpath_Computation.ipynb

# Import Libraries

# System

import os
import json

# Analysis

import numpy as np # Numerical Analysis
import geopandas as gpd
import rasterio # Rasters
from skimage import graph # least cost path analysis <- crashes due to large size

# Load universal data

cwd = os.getcwd()

raster = np.load(os.path.join(cwd,'template.npy'))

trip = gpd.read_file(os.path.join(cwd, 'DoryTrip.geojson')).to_crs('EPSG:26915')

# Definitions

def Save_Geotiff_to_template(array, template_path, save_path):
    '''Saves a numpy array into a geotiff with the same CRS as the template.
    '''            

    # Get metadata from template
    rst = rasterio.open(template_path) # Open template
    meta = rst.meta.copy() # Copy template metadata
    # meta.update(compress='lzw') # Good for integers/categorical rasters
    rst.close()

    with rasterio.open(save_path, 'w+', **meta) as out: # Burn features into raster
        out.write_band(1, array)
        
def get_lowest_cost_path(cost_surface, trip_gdf, numpy_raster_template):
    ''' Creates a raster representing the lowest cost path.
    Returns total cost of path and raster
    '''

    ## Find raster index for trip
    
    raster = numpy_raster_template

    # Start

    start_pt = trip_gdf[trip_gdf['type'] == 'start'].geometry[0]
    start_coords = [start_pt.x, start_pt.y]

    start_x_index = np.unravel_index(np.abs(raster[0] - start_coords[0]).argmin(axis = None), raster[0].shape)[0]
    start_y_index = np.unravel_index(np.abs(raster[1] - start_coords[1]).argmin(axis = None), raster[1].shape)[1]
    start_index = [start_x_index, start_y_index]

    # End

    end_pt = trip_gdf[trip_gdf['type'] == 'end'].geometry[1]
    end_coords = (end_pt.x, end_pt.y)

    end_x_index = np.unravel_index(np.abs(raster[0] - end_coords[0]).argmin(axis = None), raster[0].shape)[0]
    end_y_index = np.unravel_index(np.abs(raster[1] - end_coords[1]).argmin(axis = None), raster[1].shape)[1]
    end_index = (end_x_index, end_y_index)
    
    ## Get path
    
    #m = graph.MCP(-1*cost_surface.T, fully_connected=False)
    
    #costs, traceback_array = m.find_costs([start_index], [end_index])
    
    #indices = m.traceback(end_index)
    #cost = costs[end_index]
    
    indices, cost = graph.route_through_array(-1*cost_surface.T, start_index, end_index,
                        fully_connected=False, geometric=False)
    
    ## Create raster of the path
    
    ids = np.stack(indices, axis=-1)
    path = np.zeros_like(cost_surface.T)
    path[ids[0], ids[1]] = 1
    
    return cost, path


# Get All Filenames

file_names = np.array(os.listdir('5_weights_tests'))

# Get least cost paths

cost_of_paths = {}
    
for rast_name in np.array(file_names):

    # Load that raster

    filepath = os.path.join('5_weights_tests', rast_name)

    rast = rasterio.open(filepath) # Open

    values = rast.read(1) # Save band

    rast.close() # Close

    # Get cost path
    cost, costpath = get_lowest_cost_path(values, trip, raster)
    
    cost_of_paths[rast_name] = cost
    
    # Save as Geotiff
        
    savepath = os.path.join('6_costpaths', f'{rast_name[:-17]}_path.tif')
        
    Save_Geotiff_to_template(costpath, 'template.tif', savepath)
    
# Save paths as Json

with open("cost_of_paths.json", "w") as outfile:
    json.dump(cost_of_paths, outfile)
