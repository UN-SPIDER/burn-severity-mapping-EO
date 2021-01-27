# UN-SPIDER Recommended Practice on Burn Severity Mapping with Satellite Data
## Overview
Wildfires can result in the loss of human life and influence different ecological processes as they partially or completely remove the vegetation layer. Therefore, it is essential to assess the severity of the impacted area. This Recommended Practice was developed to help contribute to the assessment of areas affected by wildfires. It combines the use of Landsat 8 or Sentinel-2 pre- and post-fire imagery, and the Normalized Burn Ration (NBR) index. The Recommended Practice was designed specifically to assess large areas.
The Recommended Practice can be followed using the following programming languages and GIS software:  
- [QGIS with Sentinel-2](https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/Step-by-Step/QGIS-sentinel2)  
- [QGIS with Landsat 8](https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/Step-by-Step/QGIS)
- [R with Landsat 8](https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/Step-By-Step/RStudio)
- **Python with Sentinel-2**
  * **[Jupyter Notebook (in this repository)](https://github.com/UN-SPIDER/burn-severity/blob/master/burn_severity.ipynb)**
  * **[Python script (in this repository)](https://github.com/UN-SPIDER/burn-severity/blob/master/burn_severity1.py)**
- [Google Earth Engine with Landsat 8 or Sentinel-2](https://code.earthengine.google.com/b455ba8cf4b5bee822bb7ff8935e6209)  

<div align:"center">
  <img src="https://github.com/UN-SPIDER/burn-severity/blob/master/Burn_severity_overview.JPG">
</div>

------
## Methodology
The aim of this step-by-step procedure is the generation of a burn severity map for the assessment of the areas affected by wildfires. The Normalized Burn Ratio (NBR) is used, as it was designed to highlight burned areas and estimate burn severity. It uses near-infrared (NIR) and shortwave-infrared (SWIR) wavelengths. Healthy vegetation before a fire has very high NIR reflectance and a low SWIR response. In contrast, recently burned areas have a low reflectance in the NIR and high reflectance in the SWIR band. More information about the NBR can be found on the [UN-SPIDER Knowledge Portal](https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/in-detail/normalized-burn-ratio). 
The NBR is calculated for images before the fire (pre-fire NBR) and for images after the fire (post-fire NBR) and the post-fire image is subtracted from the pre-fire image to create the differenced (or delta) NBR (dNBR) image. The dNBR can be used for burn severity assessment, as areas with higher dNBR values indicate more severe damage whereas areas with negative dNBR values might show increased vegetation productivity. The dNBR values can be classified according to burn severity ranges proposed by the United States Geological Survey (USGS).  
## Data and software required
This Recommended Practice requires either Landsat 8 or Sentinel-2 optical satellite imagery. Landsat 8 data can be downloaded through the [USGS Earth Explorer](https://un-spider.org/node/10960) platform, while Sentinel-2 data is available through the [Copernicus Open Access Hub](https://un-spider.org/fr/links-and-resources/data-sources/batch-download-sentinel).
The following images are needed:
- Landsat 8/Sentinel-2 pre-fire NIR and SWIR2/SWIR images
- Landsat 8/Sentinel-2 post-fire NIR and SWIR2/SWIR images
- Shapefile of the area of interest
- QGIS, R, Python or Google Earth Engine  

The following parts will outline the approach with Python and Anaconda.

## Approach with Python
### Download and install Anaconda 3
To run the following Jupyter Notebook you need python installed on your computer. We recommend downloading and installing Anaconda 3 (https://www.anaconda.com/download/) with the python 3.6 version as it includes a lot of useful packages.

#### Preparing the working environment

There are two possible approaches:  
1. Install the required packages directly
2. Use the provided environment file to setup a new environment

##### Installing the required packages directly

The required package versions are:

- python 3.6 or 3.7
- gdal >= 2.3.1
- numpy >= 1.15.0

The gdal package needs to be installed and it must have Bigtiff support; and numpy may need to be updated. For these reasons it is recommended to use the [conda-forge package repository](https://conda-forge.org/) as follows:

1. Open the Anaconda Prompt (for instance by searching for it in the Windows start menu).

2. Type `conda install -c conda-forge --no-pin gdal numpy` and hit `Enter`.

3. Type `y` and hit `Enter` if you are asked whether you want to proceed.

**Note:** `--no-pin` is added to the install command in order to not restrict maximum package versions, and thus to avoid having certain packages be downgraded when installing gdal. If gdal and/or numpy do not meet the version requirements, then modify the install command to:

`conda install -c conda-forge --no-pin gdal=2.3.1 numpy=1.15.0`

##### Set up a new environment

All the relevant documentation can be found in the [Anaconda documentation](https://conda.io/docs/user-guide/tasks/manage-environments.html).

Carry out the following steps in order to setup the environment:

1. Download the `environment.yml` file from this repository.
2. Open the Anaconda Prompt (for instance by searching for it in the Windows start menu).
3. Navigate to where you have saved the environment file: `cd C:\my\path`
4. Type `conda env create -f environment.yml` and hit `Enter`.

After Step 1 verify that the file you have downloaded has the following structure:

**name:** burn-severity
**channels:**
conda-forge
defaults
**dependencies:**
python=3.6
gdal=2.2.4
numpy=1.15.*
matplotlib
geopandas= 0.4.0
rasterio= 1.0.7
ipython
jupyter

The ***burn-severity*** environment should now be created. Make sure it is activated before running Jupyter or before running the script if you save it as a file (see the general documentation for how to do this).

## Processing steps
In order to generate the burn severity map with Sentinel-2 imagery, the Sentinel-2 bands B8A and B12 are used. Thus, in the Python code (after the function definitions) the bands of both pre- and post-fire images are loaded and the NBR is calculated. As it is already mentioned NBR is calculated using the NIR and SWIR bands of Sentinel-2 in this practice, using the formula shown below:

`NBR = (NIR - SWIR) / (NIR + SWIR)`

`NBR = (B8A - B12) / (B8A + B12)`

After calculating NBR for both pre and post-fire images, it is possible to determine their difference and obtain dNBR. Therefore, the difference between the NBR before and after the fire, referred to as dNBR during this practice is calculated. For that, the NBR of the post-images is subtracted from the NBR of the pre-images, as shown in the formula below:

`dNBR = pre-fire NBR – post-fire NBR`

In addition, the shapefile of the correspondent area is used to clip the raster dNBR to obtain an image reflecting only the study area. As during this recommended practice, we are using the area of Empedrado region in Chile, the administrative boundary of the aforementioned area is used. In order to be able to clip the data, the shapefile and dNBR should be in the same projection and thus the shapefile is firt reprojected.

Then, USGS standard is used to classify the burn severity map according to the proposed number ranges, as explained here.

Finally, in order to quantify the area in each burn severity class, the statistics of the raster must be calculated. A transformation of the raster, so that all pixels are assigned one value for each burn severity class is necessary before the calculation.
