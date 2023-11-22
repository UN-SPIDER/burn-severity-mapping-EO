import zipfile
import geopandas as gpd
import os
import tempfile

def ingest_BARC_zip_file(zip_file_path):
    valid_shapefiles = []
    
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        for file_name in zip_ref.namelist():
            if file_name.endswith('.shp'):
                # Create a temporary directory
                with tempfile.TemporaryDirectory() as tmp_dir:
                    # Extract the related files to the temporary directory
                    shp_base = os.path.splitext(file_name)[0]
                    for ext in ['.shp', '.shx', '.dbf']:
                        zip_ref.extract(shp_base + ext, path=tmp_dir)
                    
                    print("Found shapefile: {}".format(file_name))

                    # Read the shapefile from the temporary directory
                    valid_shapefile = gpd.read_file(os.path.join(tmp_dir, file_name))
                    valid_shapefiles.append(valid_shapefile)

    return valid_shapefiles