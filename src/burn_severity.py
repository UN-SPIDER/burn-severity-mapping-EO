import os
from osgeo import osr
from osgeo import ogr
from osgeo import gdal
import numpy as np
import boto3
from botocore.exceptions import NoCredentialsError
from botocore.handlers import disable_signing
import requests
import matplotlib
import matplotlib.pyplot as plt
import rasterio
from rasterio.merge import merge
import glob
from rasterio.plot import show
from rasterio.mask import mask
from shapely.geometry import mapping
import geopandas as gpd
import math
import rioxarray as rxr
import pandas as pd

os.environ["AWS_NO_SIGN_REQUEST"] = "YES" # to avoid signing requests, avoid AWS auth

def read_band_image_from_stac_item_collection(band, item_collection):
    """
    This function takes as input the Sentinel-2 band name and a PySTAC ItemCollection,
    reads the image from each item's asset for the given band, and returns the merged
    data from all images as an array.
    input:   band           string            Sentinel-2 band name
             item_collection ItemCollection   PySTAC ItemCollection
    output:  data           array (n x m)     array of the band image
             spatialRef     string            projection 
             geoTransform   tuple             affine transformation coefficients
             targetprj                        spatial reference
    """
    images = []
    for item in item_collection:
        asset = item.assets.get(band)
        if asset is not None:
            image = rasterio.open(asset.href)
            images.append(image)
    print('Merging - number of images: ', len(images))
    mosaic, __out_trans = merge(images)
    spatialRef = images[0].crs
    geoTransform = images[0].get_transform()
    targetprj = osr.SpatialReference(wkt = images[0].crs.wkt)
    return mosaic, spatialRef, geoTransform, targetprj

def calc_nbr(band_nir, band_swir):
    """
    This function takes an input the arrays of the bands from the read_band_image
    function and returns the Normalized Burn ratio (NBR)
    input:  band_nir   array (n x m)      array of first band image e.g B8A
            band_swir   array (n x m)      array of second band image e.g. B12
    output: nbr     array (n x m)      normalized burn ratio
    """
    nbr = (band_nir - band_swir) / (band_nir + band_swir)
    return nbr

def calc_dnbr(nbr_prefire, nbr_postfire):
    """
    This function takes as input the pre- and post-fire NBR and returns the dNBR
    input:  nbr_prefire     array (n x m)       pre-fire NBR
            nbr_postfire     array (n x m)       post-fire NBR
    output: dnbr     array (n x m)       dNBR
    """
    dnbr = nbr_prefire - nbr_postfire
    return dnbr

def calc_rdnbr(dnbr, nbr_prefire):
    """
    This function takes as input the dNBR and prefire NBR, and returns the relative dNBR
    input:  dnbr     array (n x m)       dNBR
            nbr_prefire     array (n x m)       pre-fire NBR
    output: rdnbr    array (n x m)       relative dNBR
    """
    rdnbr = dnbr / np.abs(np.square(nbr_prefire))
    return rdnbr

def calc_rbr(dnbr, nbr_prefire):
    """
    This function takes as input the dNBR and prefire NBR, and returns the relative burn ratio
    input:  dnbr     array (n x m)       dNBR
            nbr_prefire     array (n x m)       pre-fire NBR
    output: rbr    array (n x m)       relative burn ratio
    """
    rbr = dnbr / (nbr_prefire + 1.001)
    return rbr

def calc_burn_metrics(prefire_nir, prefire_swir, postfire_nir, postfire_swir):
    """
    This function takes as input the pre- and post-fire NIR and SWIR bands and returns the
    NBR, dNBR, rDNBR, and rBR
    input:  prefire_nir     array (n x m)       pre-fire NIR
            prefire_swir     array (n x m)       pre-fire SWIR
            postfire_nir     array (n x m)       post-fire NIR
            postfire_swir     array (n x m)       post-fire SWIR
    output: nbr_prefire     array (n x m)       normalized burn ratio
            nbr_postfire     array (n x m)       normalized burn ratio
            dnbr     array (n x m)       dNBR
            rdnbr    array (n x m)       relative dNBR
            rbr    array (n x m)       relative burn ratio
    """
    nbr_prefire = calc_nbr(prefire_nir, prefire_swir)
    nbr_postfire = calc_nbr(postfire_nir, postfire_swir)
    dnbr = calc_dnbr(nbr_prefire, nbr_postfire)
    rdnbr = calc_rdnbr(dnbr, nbr_prefire)
    rbr = calc_rbr(dnbr, nbr_prefire)
    # stack these arrays together, naming them by their source

    burn_stack = rxr.concat(
        [nbr_prefire, nbr_postfire, dnbr, rdnbr, rbr],
        pd.Index(['nbr_prefire', 'nbr_postfire', 'dnbr', 'rdnbr', 'rbr'], name='burn_metric')
    )

    return burn_stack

def reproject_shp_gdal(infile, outfile, targetprj):
    """
    This function takes as input the input and output file names and the projection
    in which the input file will be reprojected and reprojects the input file using
    gdal
    input:  infile     string      input filename
            outfile    string      output filename
            targetprj              projection (output of function read_band_image)
    """
    ## reprojection with gdal 
    
    driver = ogr.GetDriverByName("ESRI Shapefile") 
    dataSource = driver.Open(infile, 1) # 0 means read-only. 1 means writeable.
    layer = dataSource.GetLayer()
    sourceprj = layer.GetSpatialRef()
    transform = osr.CoordinateTransformation(sourceprj, targetprj)
    
    # Create the output shapefile
    outDriver = ogr.GetDriverByName("Esri Shapefile")
    outDataSource = outDriver.CreateDataSource(outfile)
    outlayer = outDataSource.CreateLayer('', targetprj, ogr.wkbPolygon)
    outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    
    #Iterate over Features
    i = 0
    for feature in layer:
        transformed = feature.GetGeometryRef()
        transformed.Transform(transform) #reproject geometry

        geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb()) # create geometry from wkb (write geometry of reprojected geometry)
        defn = outlayer.GetLayerDefn() #layer definition
        feat = ogr.Feature(defn)  #create new feature
        feat.SetField('id', i) #set id
        feat.SetGeometry(geom) #set geometry
        outlayer.CreateFeature(feat) 
        i += 1
        feat = None
        
def array2raster(array, geoTransform, projection, filename):
    """ 
    This function tarnsforms a numpy array to a geotiff projected raster
    input:  array                       array (n x m)   input array
            geoTransform                tuple           affine transformation coefficients
            projection                  string          projection
            filename                    string          output filename
    output: dataset                                     gdal raster dataset
            dataset.GetRasterBand(1)                    band object of dataset
    
    """
    pixels_x = array.shape[1]
    pixels_y = array.shape[0]
    
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(
        filename,
        pixels_x,
        pixels_y,
        1,
        gdal.GDT_Float64, )
    dataset.SetGeoTransform(geoTransform)
    dataset.SetProjection(projection)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.
    return dataset, dataset.GetRasterBand(1)  #If you need to return, remenber to return  also the dataset because the band don`t live without dataset.
 
def clip_raster(filename, shp):
    """
    This function clips a raster based on a shapefile
    input:  filename          string                input raster filename
            shp               dataframe             input shapefile open with geopandas
    output: clipped           array (1 x n x m)     clipped array 
            clipped_meta      dict                  metadata
            cr_ext            tuple                 extent of clipped data
            gt                tuple                 affine transformation coefficients
    """
    inraster = rasterio.open(filename)
    
    extent_geojson = mapping(shp['geometry'][0])
    clipped, crop_affine = mask(inraster, 
                                shapes=[extent_geojson], 
                                nodata = np.nan,
                                crop=True)
    clipped_meta = inraster.meta.copy()
    # Update the metadata to have the new shape (x and y and affine information)
    clipped_meta.update({"driver": "GTiff",
                 "height": clipped.shape[0],
                 "width": clipped.shape[1],
                 "transform": crop_affine})
    cr_ext = rasterio.transform.array_bounds(clipped_meta['height'], 
                                            clipped_meta['width'], 
                                            clipped_meta['transform'])
    
    # transform to gdal
    gt = crop_affine.to_gdal()
    
    return clipped, clipped_meta, cr_ext, gt
    
def reclassify(array):
    """
    This function reclassifies an array
    input:  array           array (n x m)    input array
    output: reclass         array (n x m)    reclassified array
    """
    reclass = np.zeros((array.shape[0],array.shape[1]))
    for i in range(0,array.shape[0]):
        for j in range(0,array.shape[1]):
            if math.isnan(array[i,j]):
                reclass[i,j] = np.nan
            elif array[i,j] < 0.1:
                reclass[i,j] = 1
            elif array[i,j] < 0.27:
                 reclass[i,j] = 2
            elif array[i,j] < 0.44:
                 reclass[i,j] = 3
            elif array[i,j] < 0.66:
                 reclass[i,j] = 4
            else:
                reclass[i,j] = 5
                
    return reclass

def is_s3_url_valid(url):
    s3 = boto3.client('s3')
    s3.meta.events.register('choose-signer.s3.*', disable_signing)

    bucket_name = url.split('/')[2]
    key = '/'.join(url.split('/')[3:])
    try:
        response = s3.list_objects_v2(Bucket=bucket_name, Prefix=key)
        for obj in response.get('Contents', []):
            if obj['Key'] == key:
                return True
        return False
    except NoCredentialsError:
        print("No AWS credentials found")
        return False
    except Exception as e:
        print(f"Invalid S3 URL: {url}. Exception: {str(e)}")
        return False