import requests
import geopandas as gpd
import rasterio.features
from sentinelhub import WmsRequest, WcsRequest, MimeType, CRS, BBox
import numpy as np

INSTANCE_ID = 'YOUR-INSTANCE-ID'

def query_sentinelhub(input_shp, resolution, time_range):

    bounds = shp.bounds.values[0]
    bbox = BBox(bbox=list(bounds), crs=CRS.WGS84)

    wms=True

    if wms: 
        wms_request = WmsRequest(layer='TRUE-COLOR-S2-L1C',
                                bbox=bbox,
                                time=time_range,
                                width=512, height=856,
                                instance_id=INSTANCE_ID)
        img = wms_request.get_data()
    else: 
        resolution = '60m' 
        wcs_request = WcsRequest(layer='TRUE-COLOR-S2-L1C',
                                bbox=bbox,
                                time=time_range,
                                resx=resolution, resy=resolution,
                                instance_id=INSTANCE_ID)
        img = wcs_request.get_data()
        return img

    #select first image in time series
    img_median = sentinel_median_over_time(img)

    #create a rasterio dataset in memory
    with rasterio.MemoryFile() as memfile:
        with memfile.open(driver='GTiff',
                        width=img_median.shape[1], #ncols for the raster
                        height=img_median.shape[0], #nrows for the raster
                        count=1, #number of bands
                        dtype=str(img_median.dtype)) as dataset:
            dataset.write(img_median[:,:,1],1) #write numpy array to raster
            return dataset

def sentinel_median_over_time(sentinel_list):
    #concatenate images along a new time axis
    timeseries_img = np.stack(data, axis=0)

    #calculate median along time axis
    median_img = np.median(timeseries_img, axis=0)

    return median_img