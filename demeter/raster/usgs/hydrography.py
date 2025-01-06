"""
Tools for fetching hydrography data from USGS:
https://www.usgs.gov/national-hydrography/access-national-hydrography-products

Example:

    raster, transform, crs = fetch_and_merge_rasters("cat.tif", "path/to/boundaries.geojson")
"""

import os
import re
from collections.abc import Iterable, Sequence
from typing import Union
from zipfile import ZipFile

import geopandas

from demeter.raster.usgs.constants import S3_BUCKET_NAME
from demeter.raster.usgs.utils import (
    download_from_s3,
    merge_and_crop_rasters,
    s3_client,
)
from demeter.raster.utils import Raster

# USGS rasters are organized by 4-digit Hydrologic Unit (HU4). To know which
# rasters to download, we need to identify which HU4 regions the input geometry
# intersects with. To do this, I extracted the HU4 layer from USGS' Watershed
# Boundary Dataset (WBD) and removed features outside the CONUS:
# https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/WBD/National/
#
# To keep this file under GitHub's 100mb size limit, I converted it to topojson
# (https://github.com/topojson/topojson) and compressed it.
# TODO: project to EPSG:5070
HU4_CONUS_PATH = os.path.join(os.path.dirname(__file__), "hu4_conus.zip")

# USGS hydrography rasters are stored in S3 under:
S3_PREFIX = "StagedProducts/Hydrography/NHDPlusHR/VPU/Current/Raster/"

# USGS raster files for CONUS use the EPSG:5070 projection. Some rasters
# specify the ESRI:102039 projection, which is essentially the same thing.
RASTER_CRS = "EPSG:5070"


def fetch_and_merge_rasters(
    raster_filename: str,
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries],
    crop: bool = True,
) -> Raster:
    """
    Fetch the given raster (e.g. "cat.tif") from USGS for the given geometries.
    If the geometries span multiple HU4 regions, fetch all the necessary
    rasters and stitch them together.

    If `crop` is True (the default), crop the output raster to the given
    geometries.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    rasters = list(fetch_rasters(raster_filename, geometries))

    return merge_and_crop_rasters(
        rasters, crop_to=geometries.to_crs(RASTER_CRS) if crop else None
    )


def fetch_rasters(
    raster_filename: str,
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries],
) -> Iterable[str]:
    """
    Fetch all the rasters with the given filename (e.g. "cat.tif") that
    intersect with the given geometries. Yield the path to each downloaded
    raster.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    if geometries.empty:
        raise ValueError("No geometries provided")

    if not raster_filename.endswith(".tif"):
        raster_filename = f"{raster_filename}.tif"

    hu4_codes = find_hu4_codes(geometries)
    archive_paths = download_raster_archives(hu4_codes)

    for archive_path in archive_paths:
        zip_archive = ZipFile(archive_path)
        raster_paths = [
            path
            for path in zip_archive.namelist()
            if os.path.basename(path) == raster_filename
        ]

        if not raster_paths:
            raise Exception(
                f"Could not find raster '{raster_filename}' in {archive_path}"
            )
        if len(raster_paths) > 1:
            raise Exception(
                f"Multiple '{raster_filename}' files found in {archive_path}"
            )

        raster_path = raster_paths[0]
        yield f"zip://{archive_path}!{raster_path}"


def find_hu4_codes(
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries]
) -> Sequence[str]:
    """
    Return the HU4 codes for the regions that intersect with the given
    geometries.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    # TODO: project both input geometries and HU4 regions to EPSG:5070 first
    # for more accurate intersection in CONUS
    projected_geometries = geometries.geometry.to_crs(epsg=4269)
    hu4_regions = geopandas.read_file(
        HU4_CONUS_PATH,
        mask=projected_geometries.union_all(),
    )

    if hu4_regions.empty:
        raise ValueError("No HU4 regions found for geometries. Are they in CONUS?")

    geometries_without_hu4_region = geometries[
        projected_geometries.disjoint(hu4_regions.union_all())
    ]
    if not geometries_without_hu4_region.empty:
        raise ValueError(
            f"Can't find HU4 region for geometries at index: {geometries_without_hu4_region.index.tolist()}. Are they in CONUS?"
        )

    return hu4_regions["huc4"].tolist()


def download_raster_archives(hu4_codes: Iterable[str]) -> Iterable[str]:
    """
    Download the raster .zip files for the given HU4 codes.
    """
    raster_keys = raster_keys_by_hu4_code()

    for hu4_code in hu4_codes:
        key = raster_keys[hu4_code]
        yield download_from_s3(key)


_raster_keys_by_hu4_code: dict[str, str] = {}  # populate this lazily


def raster_keys_by_hu4_code() -> dict[str, str]:
    """
    Return a dict mapping HU4 codes to raster keys in S3.
    """
    if not _raster_keys_by_hu4_code:
        raster_keys = _fetch_raster_keys()

        for key in raster_keys:
            match = re.search(r"NHDPLUS_H_(\d{4})_HU4(?:_\d{8})?_RASTER\.zip$", key)
            if not match:
                continue

            hu4_code = match.group(1)
            _raster_keys_by_hu4_code[hu4_code] = key

    return _raster_keys_by_hu4_code


def _fetch_raster_keys() -> Sequence[str]:
    print("Fetching list of available rasters")
    response = s3_client.list_objects_v2(
        Bucket=S3_BUCKET_NAME,
        Prefix=S3_PREFIX,
    )

    return [
        item["Key"] for item in response["Contents"] if item["Key"].endswith(".zip")
    ]
