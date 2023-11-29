import requests
import geopandas as gpd
from pystac_client import Client
from datetime import datetime
import planetary_computer

SENTINEL2_PATH = "https://planetarycomputer.microsoft.com/api/stac/v1"


class Sentinel2Client:
    def __init__(self, geojson_bounds, buffer = .1):
        self.path = SENTINEL2_PATH
        self.client = Client.open(
            self.path,
            modifiers = planetary_computer.sign_inplace
        )
        self.geojson_bounds = geojson_bounds
        geojson_bbox = geojson_bounds.bounds.to_numpy()[0]
        self.bbox = [
            geojson_bbox[0].round(decimals=2) - buffer,
            geojson_bbox[1].round(decimals=2) - buffer,
            geojson_bbox[2].round(decimals=2) + buffer,
            geojson_bbox[3].round(decimals=2) + buffer
        ]



    def query(self, date_range, cloud_cover=20, from_bbox = False, max_items = None):

        date_range_fmt = "{}/{}".format(
            date_range[0], date_range[1]
        )

        query = {
            "collections": ["sentinel-2-l2a"],
            "datetime": date_range_fmt,
            # "query": {"eo:cloud_cover": {"lt": cloud_cover}}
        }

        if from_bbox:
            query["bbox"] = self.bbox
        else:
            query["intersects"] = self.geojson_bounds

        if max_items:
            query["max_items"] = max_items

        items = self.client.search(**query).item_collection()

        return items
