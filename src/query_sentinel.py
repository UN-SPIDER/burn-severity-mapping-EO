import requests
import geopandas as gpd
from pystac_client import Client
from datetime import datetime

SENTINEL2_PATH = "https://earth-search.aws.element84.com/v1"


class Sentinel2Client:
    def __init__(self, geojson_bounds):
        self.path = SENTINEL2_PATH
        self.client = Client.open(self.path)
        self.geojson_bounds = geojson_bounds

    def query(self, date_range, cloud_cover=20):

        date_range_fmt = "{}/{}".format(
            date_range[0], date_range[1]
        )

        query = {
            "collections": ["sentinel-s2-l2a-cogs"],
            "intersects": self.geojson_bounds,
            "datetime": date_range_fmt,
            "query": {"eo:cloud_cover": {"lt": cloud_cover}},
        }

        items = self.client.search(**query).get_all_items()

        return items
