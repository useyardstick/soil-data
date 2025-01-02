import os

# This is where USGS keeps its raster archives:
S3_BUCKET_NAME = "prd-tnm"

# Downloaded raster files are cached here:
CACHED_RASTER_FILES_DIRECTORY = os.environ.get(
    "USGS_CACHED_RASTER_FILES_DIRECTORY", ".usgs_cache"
)
