import calendar
from collections.abc import Iterable, Sequence

import requests

from demeter.raster.sentinel2.constants import ODATA_PRODUCTS_ENDPOINT, S3_BUCKET_NAME
from demeter.raster.sentinel2.utils.tiles import TileMetadata


def find_safe_files(
    tiles: Iterable[TileMetadata],
    year: int,
    month: int,
) -> Iterable[str]:
    """
    Yield the S3 keys for Sentinel-2 SAFE files for the given (tile_id,
    orbit_number) pairs recorded during the given month.
    """
    for tile in tiles:
        yield from find_safe_files_for_tile(tile, year, month)


def find_safe_files_for_tile(
    tile: TileMetadata,
    year: int,
    month: int,
) -> Sequence[str]:
    """
    Using the Copernicus OData search endpoint, find all Sentinel-2 SAFE files
    for the given (tile_id, orbit_number) pair recorded during the given
    month, and return their S3 keys. This is much faster than scanning the
    bucket.

    See https://documentation.dataspace.copernicus.eu/APIs/OData.html
    """
    query = _odata_query(tile, year, month)
    limit = 100
    response = requests.get(
        ODATA_PRODUCTS_ENDPOINT, params={"$filter": query, "$top": str(limit)}
    )
    response.raise_for_status()
    s3_keys = [
        item["S3Path"].removeprefix(f"/{S3_BUCKET_NAME}/")
        for item in response.json()["value"]
    ]

    # There are typically only 6 results per month for any given (tile, orbit)
    # pair, so we shouldn't ever need to paginate these results. Just to be
    # safe, make sure we didn't hit the result limit:
    assert len(s3_keys) < limit
    assert all(key.endswith(".SAFE") for key in s3_keys)
    return s3_keys


def _odata_query(
    tile: TileMetadata,
    year: int,
    month: int,
) -> str:
    assert 1 <= month <= 12
    tile_id, relative_orbit_number = tile
    _, last_day_of_month = calendar.monthrange(year, month)
    conditions = [
        "Collection/Name eq 'SENTINEL-2'",
        f"ContentDate/Start ge {year}-{month:02}-01",
        f"ContentDate/End le {year}-{month:02}-{last_day_of_month:02}",
        "Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq 'S2MSI2A')",
        f"Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'tileId' and att/OData.CSC.StringAttribute/Value eq '{tile_id}')",
        f"Attributes/OData.CSC.IntegerAttribute/any(att:att/Name eq 'relativeOrbitNumber' and att/OData.CSC.IntegerAttribute/Value eq {relative_orbit_number})",
    ]
    return " and ".join(conditions)
