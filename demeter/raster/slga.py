"""
Tools for fetching raster data from the Soil and Landscape Grid of Australia
(SLGA): https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html

Example:

    slga_om = fetch_slga_data_for_depth_range(
        geometries=[
            {
                "type": "Polygon",
                "coordinates": [...],
            },
            ...
        ],
        soil_property=SoilProperty.ORGANIC_CARBON,
        start_depth=0,
        end_depth=100,
    )
    mean_raster, transform, crs = slga_om.mean
    stddev_raster, _, _ = slga_om.stddev
"""

__all__ = [
    "Depth",
    "SoilProperty",
    "Statistic",
    "estimate_carbon_stock",
    "fetch_slga_data_for_depth_range",
]

import itertools
import math
import os
from dataclasses import dataclass
from enum import Enum, unique
from typing import Optional, Union
from urllib.parse import urljoin

import geopandas
import rasterio
import rasterio.transform
import rasterio.windows

from demeter.raster.utils import (
    DepthEnum,
    Raster,
    calculate_carbon_stock_stddev,
    calculate_weighted_average_mean,
    calculate_weighted_average_stddev,
    mask_raster,
)

BASE_URL = "https://apikey:{tern_api_key}@data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/"


@dataclass(frozen=True)
class CombinedRasters:
    mean: Raster
    stddev: Optional[Raster] = None


@unique
class SoilProperty(Enum):
    BULK_DENSITY = "BDW/BDW_{start_depth:03d}_{end_depth:03d}_{statistic}_N_P_AU_TRN_N_20230607.tif"
    ORGANIC_CARBON = "SOC/SOC_{start_depth:03d}_{end_depth:03d}_{statistic}_N_P_AU_TRN_N_20220727.tif"
    CLAY = "CLY/CLY_{start_depth:03d}_{end_depth:03d}_{statistic}_N_P_AU_TRN_N_20210902.tif"
    PH = "PHW/PHW_{start_depth:03d}_{end_depth:03d}_{statistic}_N_P_AU_TRN_N_20220520.tif"


@unique
class Statistic(Enum):
    MEAN = "EV"
    P5 = "05"
    P95 = "95"


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
    bulk density rasters from SLGA, then combining them into a estimated carbon
    stock raster.
    """
    organic_carbon, bulk_density = [
        fetch_slga_data_for_depth_range(
            geometries,
            soil_property=soil_property,
            start_depth=start_depth,
            end_depth=end_depth,
            calculate_standard_deviation=calculate_standard_deviation,
        )
        for soil_property in [SoilProperty.ORGANIC_CARBON, SoilProperty.BULK_DENSITY]
    ]

    assert organic_carbon.mean.shape == bulk_density.mean.shape

    # SLGA organic carbon and bulk density GeoTIFF files don't have exactly
    # the same dimensions. They should still be very close though, so it
    # should be fine to use the same transform for both.
    assert organic_carbon.mean.transform.almost_equals(bulk_density.mean.transform)
    transform = organic_carbon.mean.transform

    soil_organic_carbon_mean = organic_carbon.mean.pixels
    bulk_density_mean = bulk_density.mean.pixels
    carbon_stock_mean = soil_organic_carbon_mean * bulk_density_mean

    if not calculate_standard_deviation:
        return CombinedRasters(mean=Raster(carbon_stock_mean, transform, "EPSG:4326"))

    assert organic_carbon.stddev
    assert bulk_density.stddev
    assert organic_carbon.stddev.shape == bulk_density.stddev.shape
    assert organic_carbon.stddev.transform.almost_equals(bulk_density.stddev.transform)
    soil_organic_carbon_stddev = organic_carbon.stddev.pixels
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


def fetch_slga_data_for_depth_range(
    geometries: Union[str, geopandas.GeoSeries, geopandas.GeoDataFrame],
    soil_property: Union[str, SoilProperty],
    *,
    start_depth: int = 0,
    end_depth: int,
    calculate_standard_deviation: bool = True,
) -> CombinedRasters:
    """
    High-level interface to SLGA.

    Fetch all SLGA pixels between the given `start_depth` and `end_depth`, and
    return a depth-weighted average across the entire depth range.

    If `calculate_standard_deviation` is True (default), also return a raster
    showing the standard deviation of the distribution at each pixel.
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
    #             Statistic.MEAN: ...,
    #             Statistic.P5: ...,
    #             Statistic.P95: ...,
    #         },
    #         Depth.FIVE_TO_FIFTEEN_CM: ...,
    #         ...
    #     }
    bounds = tuple(geometries.total_bounds.tolist())
    depths = Depth.select_between(start_depth, end_depth)
    statistics = [Statistic.MEAN]
    if calculate_standard_deviation:
        statistics += [Statistic.P5, Statistic.P95]

    rasters_by_depth: dict[Depth, dict[Statistic, Raster]] = {}
    for depth, statistic in itertools.product(depths, statistics):
        if depth not in rasters_by_depth:
            rasters_by_depth[depth] = {}
        rasters_by_depth[depth][statistic] = _download_slga_raster(
            bounds, soil_property, statistic, depth
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

    # Calculate weighted average for depth range:
    weights = [depth.thickness for depth in depths]
    mean_rasters = [rasters_by_depth[depth][Statistic.MEAN].pixels for depth in depths]
    weighted_average_mean = calculate_weighted_average_mean(mean_rasters, weights)

    # Crop to geometries:
    cropped_mean_raster = mask_raster(
        Raster(weighted_average_mean, transform, "EPSG:4326"),
        geometries,
        crop=True,
        all_touched=True,
    )

    if not calculate_standard_deviation:
        return CombinedRasters(mean=cropped_mean_raster)

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
        crop=True,
        all_touched=True,
    )

    return CombinedRasters(
        mean=cropped_mean_raster,
        stddev=cropped_stddev_raster,
    )


def _download_slga_raster(
    bounds: tuple[float, float, float, float],
    soil_property: SoilProperty,
    statistic: Statistic,
    depth: Depth,
    save_test_fixture: bool = False,
) -> Raster:
    """
    Download pixels for the given bounds in parallel, using the global thread pool.
    """
    path = _slga_path(soil_property, statistic, depth)

    base_url = BASE_URL.format(tern_api_key=os.environ["TERN_API_KEY"])
    url = urljoin(base_url, path)
    left, bottom, right, top = bounds

    sanitized_url = urljoin(BASE_URL, path)
    print(f"Downloading from {sanitized_url}")

    with rasterio.open(url) as dataset:
        miny, minx = rasterio.transform.rowcol(
            dataset.transform, left, top, op=math.floor
        )
        maxy, maxx = rasterio.transform.rowcol(
            dataset.transform, right, bottom, op=math.ceil
        )
        window = rasterio.windows.Window.from_slices((miny, maxy), (minx, maxx))
        pixels = dataset.read(1, window=window, masked=True)
        transform = dataset.window_transform(window)
        raster = Raster(pixels, transform, "EPSG:4326")

    # To avoid downloading from the real SLGA dataset in tests, we use local
    # test fixtures. To download a fixutre for a test, set `save_test_fixture`
    # to True and run the test once using a real Tern API key. Then set
    # `save_test_fixture` back to False, add the `mock_slga` fixture to your
    # test, and you should be good to go.
    if save_test_fixture:
        fixture_path = os.path.join(
            "tests/raster/fixtures/slga", os.path.basename(path)
        )
        raster.save(fixture_path, masked=False)

    return raster


def _slga_path(
    soil_property: SoilProperty,
    statistic: Statistic,
    depth: Depth,
) -> str:
    # HACK: One of the PH raster filenames doesn't follow the naming convention.
    # https://data.tern.org.au/landscapes/slga/NationalMaps/SoilAndLandscapeGrid/PHW/
    if (
        soil_property == SoilProperty.PH
        and statistic == Statistic.MEAN
        and depth == Depth.ZERO_TO_FIVE_CM
    ):
        return "PHW/PHW_000_005_EV_N_P_AU_TRN_N_20222005.tif"

    start_depth, end_depth = depth.value
    return soil_property.value.format(
        start_depth=start_depth,
        end_depth=end_depth,
        statistic=statistic.value,
    )
