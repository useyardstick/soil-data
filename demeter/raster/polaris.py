"""
Tools for fetching raster data from the POLARIS dataset:
http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/

Example:

```python
polaris_om = fetch_polaris_data_for_depth_range(
    "/path/to/geometries.geojson",
    soil_property="om",
    start_depth=0,
    end_depth=100,
)
mean_raster, transform, crs = polaris_om.mean
stddev_raster, _, _ = polaris_om.stddev
```
"""

__all__ = [
    "Depth",
    "SoilProperty",
    "Statistic",
    "estimate_carbon_stock",
    "fetch_polaris_data",
    "fetch_polaris_data_for_depth_range",
]

import itertools
import os
import shutil
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum, unique
from typing import Optional, Union

import geopandas
import numpy
import smart_open

from demeter.constants import OM_TO_SOC
from demeter.raster import Raster
from demeter.raster.depth_enum import DepthEnum
from demeter.raster.utils.mask import mask_raster
from demeter.raster.utils.merge import merge
from demeter.utils import (
    bounds_snapped_to_grid,
    calculate_carbon_stock_stddev,
    calculate_weighted_average_mean,
    calculate_weighted_average_stddev,
)

BASE_URL = os.environ.get(
    "POLARIS_BASE_URL", "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/"
)

# POLARIS tiles are 1° by 1°
PolarisTile = tuple[int, int, int, int]


@dataclass(frozen=True)
class CombinedRasters:
    mean: Raster
    stddev: Optional[Raster] = None
    median: Optional[Raster] = None
    mode: Optional[Raster] = None
    p5: Optional[Raster] = None
    p95: Optional[Raster] = None


@unique
class SoilProperty(Enum):
    BULK_DENSITY = "bd"
    ORGANIC_MATTER = "om"
    CLAY = "clay"
    SAND = "sand"
    SILT = "silt"
    PH = "ph"
    THETA_S = "theta_s"
    THETA_R = "theta_r"
    KSAT = "ksat"
    LAMBDA = "lambda"
    HB = "hb"
    N = "n"
    ALPHA = "alpha"


@unique
class Statistic(Enum):
    MEAN = "mean"
    MEDIAN = "p50"
    P5 = "p5"
    P95 = "p95"
    MODE = "mode"


@unique
class Depth(DepthEnum):
    ZERO_TO_FIVE_CM = (0, 5)
    FIVE_TO_FIFTEEN_CM = (5, 15)
    FIFTEEN_TO_THIRTY_CM = (15, 30)
    THIRTY_TO_SIXTY_CM = (30, 60)
    SIXTY_TO_ONE_HUNDRED_CM = (60, 100)
    ONE_HUNDRED_TO_TWO_HUNDRED_CM = (100, 200)


def estimate_carbon_stock(
    geometries: Union[str, geopandas.GeoSeries, geopandas.GeoDataFrame],
    *,
    start_depth: int = 0,
    end_depth: int,
    calculate_standard_deviation: bool = True,
) -> CombinedRasters:
    """
    Convenience function for the common use case of fetching organic matter and
    bulk density rasters from POLARIS, then combining them into a estimated
    carbon stock raster.
    """
    organic_matter, bulk_density = (
        fetch_polaris_data_for_depth_range(
            geometries,
            soil_property=soil_property,
            start_depth=start_depth,
            end_depth=end_depth,
            calculate_standard_deviation=calculate_standard_deviation,
        )
        for soil_property in [SoilProperty.ORGANIC_MATTER, SoilProperty.BULK_DENSITY]
    )

    assert organic_matter.mean.transform == bulk_density.mean.transform
    transform = organic_matter.mean.transform

    soil_organic_carbon_mean = organic_matter.mean.pixels * OM_TO_SOC
    bulk_density_mean = bulk_density.mean.pixels
    carbon_stock_mean = soil_organic_carbon_mean * bulk_density_mean

    if not calculate_standard_deviation:
        return CombinedRasters(mean=Raster(carbon_stock_mean, transform, "EPSG:4326"))

    assert organic_matter.stddev
    assert bulk_density.stddev
    assert organic_matter.stddev.transform == bulk_density.stddev.transform == transform
    soil_organic_carbon_stddev = organic_matter.stddev.pixels * OM_TO_SOC
    bulk_density_stddev = bulk_density.stddev.pixels
    carbon_stock_stddev = calculate_carbon_stock_stddev(
        soil_organic_carbon_mean,
        soil_organic_carbon_stddev,
        bulk_density_mean,
        bulk_density_stddev,
    )

    return CombinedRasters(
        mean=Raster(carbon_stock_mean, transform, "EPSG:4326"),
        stddev=Raster(carbon_stock_stddev, transform, "EPSG:4326"),
    )


def fetch_polaris_data_for_depth_range(
    geometries: Union[str, geopandas.GeoSeries, geopandas.GeoDataFrame],
    soil_property: Union[str, SoilProperty],
    *,
    start_depth: int = 0,
    end_depth: int,
    calculate_standard_deviation: bool = True,
    additional_statistics: list[Statistic] = [],
) -> CombinedRasters:
    """
    High-level interface to POLARIS.

    Fetch all POLARIS tiles between the given `start_depth` and `end_depth`,
    and return a depth-weighted average across the entire depth range.

    If `calculate_standard_deviation` is True (default), also return a raster
    showing the standard deviation at each pixel, inferred from the p5-p95
    split (assuming normal distribution).
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    if len(geometries) == 0:
        raise ValueError("Must provide at least one geometry")

    soil_property = SoilProperty(soil_property)

    # Download mean, p5 and p95 tiles for all depth layers down to the given
    # `end_depth`. Organize them in a dictionary by depth and statistic, e.g.:
    #
    #     {
    #         Depth.ZERO_TO_FIVE_CM: {
    #             Statistic.MEAN: [
    #                 "/path/to/tile_1.tif",
    #                 "/path/to/tile_2.tif",
    #                 ...,
    #             ],
    #             Statistic.P5: [...],
    #             Statistic.P95: [...],
    #         },
    #         Depth.FIVE_TO_FIFTEEN_CM: ...,
    #         ...
    #     }
    tiles = _polaris_tiles_for_geometries(geometries)
    depths = Depth.select_including(start_depth, end_depth)
    statistics = {Statistic.MEAN}
    statistics.update(additional_statistics)

    if calculate_standard_deviation:
        statistics.update([Statistic.P5, Statistic.P95])

    downloads_by_depth: dict[Depth, dict[Statistic, list[str]]] = {}
    for tile, depth, statistic in itertools.product(tiles, depths, statistics):
        if depth not in downloads_by_depth:
            downloads_by_depth[depth] = {}
        if statistic not in downloads_by_depth[depth]:
            downloads_by_depth[depth][statistic] = []
        downloads_by_depth[depth][statistic].append(
            _download_polaris_tile(tile, soil_property, statistic, depth)
        )

    # Merge tiles:
    bounds = tuple(geometries.total_bounds.tolist())
    rasters_by_depth: dict[Depth, dict[Statistic, Raster]] = {}
    for depth, downloads_by_statistic in downloads_by_depth.items():
        rasters_by_depth[depth] = {}
        for statistic, downloads in downloads_by_statistic.items():
            rasters_by_depth[depth][statistic] = merge(
                downloads, bounds=bounds, target_aligned_pixels=True
            )

    # Assert rasters all have the same shape and affine transform:
    unique_shapes = set()
    unique_transforms = set()
    for rasters_by_statistic in rasters_by_depth.values():
        for raster in rasters_by_statistic.values():
            unique_shapes.add(raster.shape)
            unique_transforms.add(raster.transform)
    assert len(unique_shapes) == len(unique_transforms) == 1
    transform = unique_transforms.pop()

    # Organic matter raster data is in log10%. Convert to decimal percentage:
    if soil_property == SoilProperty.ORGANIC_MATTER:
        for rasters in rasters_by_depth.values():
            for statistic, raster in rasters.items():
                rasters[statistic] = Raster((10**raster.pixels), transform, "EPSG:4326")

    # POLARIS sometimes returns bulk density values of -9999. Remove these:
    if soil_property == SoilProperty.BULK_DENSITY:
        for rasters in rasters_by_depth.values():
            for raster in rasters.values():
                raster.pixels[raster.pixels < 0] = numpy.ma.masked

    # Calculate weighted average for depth range:
    weights = [depth.thickness for depth in depths]

    # Adjust weights if selected depth range doesn't correspond exactly to
    # POLARIS depths:
    if start_depth > depths[0].start_depth:
        weights[0] = depths[0].end_depth - start_depth
    if end_depth < depths[-1].end_depth:
        weights[-1] = end_depth - depths[-1].start_depth

    mean_rasters = [rasters_by_depth[depth][Statistic.MEAN].pixels for depth in depths]
    weighted_average_mean = calculate_weighted_average_mean(mean_rasters, weights)

    # Crop to geometries:
    cropped_mean_raster = mask_raster(
        Raster(weighted_average_mean, transform, "EPSG:4326"),
        geometries,
        all_touched=True,
    )

    # Calculate weighted average for depth range for additional statistics:
    cropped_statistic_rasters = {}
    if len(additional_statistics) > 0:
        for statistic in additional_statistics:
            statistic_rasters = [
                rasters_by_depth[depth][statistic].pixels for depth in depths
            ]
            weighted_average_statistic = calculate_weighted_average_mean(
                statistic_rasters, weights
            )
            cropped_statistic_raster = mask_raster(
                Raster(weighted_average_statistic, transform, "EPSG:4326"),
                geometries,
                all_touched=True,
            )
            cropped_statistic_rasters[statistic] = cropped_statistic_raster

    if not calculate_standard_deviation:
        return CombinedRasters(
            mean=cropped_mean_raster,
            median=cropped_statistic_rasters.get(Statistic.MEDIAN, None),
            mode=cropped_statistic_rasters.get(Statistic.MODE, None),
            p5=cropped_statistic_rasters.get(Statistic.P5, None),
            p95=cropped_statistic_rasters.get(Statistic.P95, None),
        )

    # Calculate variance by depth using P5 and P95:
    p5_rasters = [rasters_by_depth[depth][Statistic.P5].pixels for depth in depths]
    p95_rasters = [rasters_by_depth[depth][Statistic.P95].pixels for depth in depths]
    weighted_average_stddev = calculate_weighted_average_stddev(
        p5_rasters, p95_rasters, weights
    )

    # Crop to geometries:
    cropped_stddev_raster = mask_raster(
        Raster(weighted_average_stddev, transform, "EPSG:4326"),
        geometries,
        all_touched=True,
    )

    return CombinedRasters(
        mean=cropped_mean_raster,
        stddev=cropped_stddev_raster,
        median=cropped_statistic_rasters.get(Statistic.MEDIAN, None),
        mode=cropped_statistic_rasters.get(Statistic.MODE, None),
        p5=cropped_statistic_rasters.get(Statistic.P5, None),
        p95=cropped_statistic_rasters.get(Statistic.P95, None),
    )


def fetch_polaris_data(
    geometries: Union[str, geopandas.GeoSeries, geopandas.GeoDataFrame],
    soil_property: Union[str, SoilProperty],
    statistic: Statistic,
    depth: Depth,
) -> Raster:
    """
    Low-level interface to POLARIS.

    Download raster images from POLARIS, merge them, and return a raster
    containing the values from the merged images.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    assert isinstance(geometries, (geopandas.GeoSeries, geopandas.GeoDataFrame))

    if len(geometries) == 0:
        raise ValueError("Must provide at least one geometry")

    soil_property = SoilProperty(soil_property)

    raster_paths = [
        _download_polaris_tile(tile, soil_property, statistic, depth)
        for tile in _polaris_tiles_for_geometries(geometries)
    ]
    bounds = tuple(geometries.total_bounds.tolist())
    raster = merge(raster_paths, bounds=bounds, target_aligned_pixels=True)
    return mask_raster(raster, geometries, all_touched=True)


def _polaris_tiles_for_geometries(
    geometries: Union[geopandas.GeoSeries, geopandas.GeoDataFrame]
) -> Iterable[PolarisTile]:
    """
    All the polaris tiles needed for the given geometries.
    """
    bounds = bounds_snapped_to_grid(geometries)
    tiles = itertools.chain.from_iterable(
        itertools.product(range(minx, maxx), range(miny, maxy))
        for minx, miny, maxx, maxy in bounds.itertuples(index=False)
    )
    for minx, miny in set(tiles):
        yield (minx, miny, minx + 1, miny + 1)


def _download_polaris_tile(
    tile: PolarisTile,
    soil_property: SoilProperty,
    statistic: Statistic,
    depth: Depth,
) -> str:
    path = _polaris_tile_path(tile, soil_property, statistic, depth)

    # First, check if we have a local copy of the tile:
    cache_directory = os.environ.get(
        "POLARIS_CACHED_RASTER_FILES_DIRECTORY", ".polaris_cache"
    )
    local_path = os.path.join(cache_directory, path)
    if os.path.exists(local_path):
        print(f"Cache hit: {local_path}")
        return local_path

    # If not, try to download from the remote cache, if configured:
    # NOTE: Partial files shouldn't show up in S3, so I don't think there's a race condition here.
    # https://stackoverflow.com/questions/38173710/amazon-s3-can-clients-see-the-file-before-upload-is-complete
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    remote_cache = os.environ.get("POLARIS_REMOTE_CACHE", None)
    remote_cache_path = (
        f"{remote_cache.removesuffix('/')}/{path}" if remote_cache else None
    )
    if remote_cache_path:
        try:
            with (
                smart_open.open(remote_cache_path, "rb") as remote_cache_file,
                open(local_path, "wb") as local_file,
            ):
                shutil.copyfileobj(remote_cache_file, local_file)
        except OSError:
            print(f"Remote cache miss: {remote_cache_path}")
        else:
            print(f"Downloaded from remote cache: {remote_cache_path}")
            return local_path

    # If we don't have a copy in our remote cache, download from source:
    remote_path = f"{BASE_URL.removesuffix('/')}/{path}"
    print(f"Downloading {remote_path}")
    with (
        smart_open.open(remote_path, "rb") as remote_file,
        open(local_path, "wb") as local_file,
    ):
        shutil.copyfileobj(remote_file, local_file)

    # Upload to the remote cache:
    if remote_cache_path:
        print(f"Uploading to remote cache: {remote_cache_path}")
        with (
            open(local_path, "rb") as local_file,
            smart_open.open(remote_cache_path, "wb") as remote_cache_file,
        ):
            shutil.copyfileobj(local_file, remote_cache_file)

    return local_path


def _polaris_tile_path(
    tile: PolarisTile,
    soil_property: SoilProperty,
    statistic: Statistic,
    depth: Depth,
) -> str:
    longitude_min, latitude_min, longitude_max, latitude_max = tile
    start_depth, end_depth = depth.value
    filename = f"lat{latitude_min}{latitude_max}_lon{longitude_min}{longitude_max}.tif"
    return (
        f"{soil_property.value}/{statistic.value}/{start_depth}_{end_depth}/{filename}"
    )
