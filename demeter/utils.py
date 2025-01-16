from typing import Union

import geopandas
import numpy
import pandas


def bounds_snapped_to_grid(
    geometries: Union[geopandas.GeoSeries, geopandas.GeoDataFrame],
    base: int = 1,
) -> pandas.DataFrame:
    """
    Return a (minx, miny, maxx, maxy) DataFrame with bounds encompassing the
    given geometries, snapped to a 1 degree x 1 degree grid.

    To snap to a grid of a larger size, pass a `base` argument to round
    to the nearest multiple of `base`.
    """
    # Explode geometries to avoid including unused areas in the event we get
    # multipart geometries that are very far apart:
    bounds = geometries.explode().bounds / base

    return (
        bounds.assign(
            minx=numpy.floor(bounds["minx"]),
            miny=numpy.floor(bounds["miny"]),
            maxx=numpy.ceil(bounds["maxx"]),
            maxy=numpy.ceil(bounds["maxy"]),
        )
        .astype(int)
        .drop_duplicates(ignore_index=True)
    ) * base


def calculate_carbon_stock_stddev(
    soil_organic_carbon_mean,
    soil_organic_carbon_stddev,
    bulk_density_mean,
    bulk_density_stddev,
):
    """
    Formula to combine variances:

        Var(XY) = Var(X) * Var(Y) + Var(X) * E(Y)² + Var(Y) * E(X)²

    Note that this assumes that bulk density and soil organic carbon are
    completely independent, which isn't strictly true.
    """
    soil_organic_carbon_variance = soil_organic_carbon_stddev**2
    bulk_density_variance = bulk_density_stddev**2
    carbon_stock_variance = (
        soil_organic_carbon_variance * bulk_density_variance
        + soil_organic_carbon_variance * bulk_density_mean**2
        + bulk_density_variance * soil_organic_carbon_mean**2
    )
    return numpy.sqrt(carbon_stock_variance)


def calculate_weighted_average_mean(mean_rasters, weights):
    stacked_mean = numpy.stack(mean_rasters)
    return numpy.average(stacked_mean, axis=0, weights=numpy.array(weights))


def calculate_weighted_average_stddev(p5_rasters, p95_rasters, weights):
    """
    Calculate standard deviation by depth using P5 and P95.
    """
    quantile = 1.645
    variance_rasters = []
    for p5, p95 in zip(p5_rasters, p95_rasters):
        standard_deviation = (p5 - p95) / (quantile * 2)
        variance_rasters.append(standard_deviation**2)

    stacked_variance = numpy.stack(variance_rasters)
    weighted_variance = numpy.average(
        stacked_variance, axis=0, weights=(numpy.array(weights) ** 2)
    )
    return numpy.sqrt(weighted_variance)
