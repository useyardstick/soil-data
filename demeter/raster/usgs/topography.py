"""
Tools for fetching elevation data from USGS at 1/3 arc-second resolution:
https://data.usgs.gov/datacatalog/data/USGS:3a81321b-c153-416f-98b7-cc8e5f0e17c3

Example:

```python
raster, transform, crs = fetch_and_merge_rasters("path/to/boundaries.geojson")
```
"""

import itertools
from typing import Iterable, Union

import geopandas

from demeter.raster import Raster
from demeter.raster.usgs.utils import download_from_s3, merge_and_crop_rasters
from demeter.raster.utils import bounds_snapped_to_grid

# USGS topography rasters at 1/3 arc-second resolution are stored in S3 under:
S3_PREFIX = "StagedProducts/Elevation/13/TIFF/current/"

# USGS topography rasters are in EPSG:4269. This is almost the same as
# EPSG:4326 (WGS84), but not quite!
RASTER_CRS = "EPSG:4269"


def fetch_and_merge_rasters(
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries],
    crop: bool = True,
) -> Raster:
    """
    Fetch 1/3 arc-second resolution elevation data for the given geometries
    from USGS. If the geometries span multiple 1 degree x 1 degree tiles, fetch
    all the necessary tiles and stitch them together.

    If `crop` is True (the default), crop the output raster to the given
    geometries.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    rasters = list(fetch_rasters(geometries))

    return merge_and_crop_rasters(
        rasters, crop_to=geometries.to_crs(RASTER_CRS) if crop else None
    )


def fetch_rasters(
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries],
) -> Iterable[str]:
    """
    Fetch all the 1 degree x 1 degree tiles from USGS that overlap with the
    given geometries. Yield the path to each downloaded raster.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    # FIXME: .to_crs(epsg=4269) doesn't seem to do anything for WGS84 inputs.
    projected_geometries = geometries.geometry.to_crs(RASTER_CRS)
    bounds = bounds_snapped_to_grid(projected_geometries)
    tiles = itertools.chain.from_iterable(
        itertools.product(range(minx, maxx), range(miny, maxy))
        for minx, miny, maxx, maxy in bounds.itertuples(index=False)
    )

    # Download and open the rasters for each tile:
    for minx, miny in set(tiles):
        raster_path = download_raster(
            northest_latitude=miny + 1,
            westest_longitude=minx,
        )
        yield raster_path


def download_raster(
    northest_latitude: int,
    westest_longitude: int,
) -> str:
    lat = (
        f"n{northest_latitude:02}"
        if northest_latitude >= 0
        else f"s{-northest_latitude:02}"
    )
    lng = (
        f"w{-westest_longitude:03}"
        if westest_longitude <= 0
        else f"e{westest_longitude:03}"
    )
    key = f"{S3_PREFIX}{lat}{lng}/USGS_13_{lat}{lng}.tif"
    return download_from_s3(key)
