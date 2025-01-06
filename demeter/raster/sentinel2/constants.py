import os
from enum import Enum, IntEnum, unique

# Sentinel-2 tiles, downloaded from https://sentiwiki.copernicus.eu/web/s2-products.
# Limited to tiles within CONUS only, to keep file size small.
TILES_CONUS_PATH = os.path.join(os.path.dirname(__file__), "tiles_conus.geojson")

# Sentinel-2 orbit tracks, downloaded from https://sentiwiki.copernicus.eu/web/s2-mission.
# Limited to orbits over CONUS only, to keep file size small.
ORBITS_CONUS_PATH = os.path.join(os.path.dirname(__file__), "orbits_conus.geojson")
SWATH_WIDTH_KM = 290

# Search API for Sentinel-2 data:
ODATA_PRODUCTS_ENDPOINT = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"

# Sentinel-2 data is available via S3 API:
# https://documentation.dataspace.copernicus.eu/APIs/S3.html
# Note that while the data can be accessed using S3 clients, it is not
# actually stored in AWS. Some things don't work exactly the same way, for
# example: when listing the contents of a bucket by prefix, AWS S3 lists *all*
# keys with that prefix, whereas this implementation only lists keys within a
# directory.
S3_BUCKET_NAME = "eodata"


@unique
class Band(Enum):
    BLUE = "B02"
    RED = "B04"
    NIR = "B08"
    SCL = "SCL"


@unique
class Resolution(Enum):
    R10 = "10m"
    R20 = "20m"
    R60 = "60m"


@unique
class SceneClassification(IntEnum):
    """
    What SCL band raster values mean.
    """

    CLOUD_SHADOWS = 3
    CLOUD_MEDIUM_PROBABILITY = 8
    CLOUD_HIGH_PROBABILITY = 9
    THIN_CIRRUS = 10


CLOUD_VALUES = [
    SceneClassification.CLOUD_SHADOWS,
    SceneClassification.CLOUD_MEDIUM_PROBABILITY,
    SceneClassification.CLOUD_HIGH_PROBABILITY,
    SceneClassification.THIN_CIRRUS,
]
