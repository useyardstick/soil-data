from collections.abc import Sequence

import numpy
from rasterio.merge import copy_count, copy_sum

from demeter.raster import Raster
from demeter.raster.utils import merge


def merge_min(rasters: Sequence[str]) -> Raster:
    return merge(rasters, method="min")


def merge_max(rasters: Sequence[str]) -> Raster:
    return merge(rasters, method="max")


def merge_mean(rasters: Sequence[str]) -> Raster:
    """
    Stack two numpy arrays, the first with the sum of valid values at each
    pixel, and the second with the count of valid values at each pixel. Then
    divide the first array by the second to get the mean.
    """
    raster = merge(
        rasters,
        method=_copy_sum_and_count,
        output_count=2,
    )
    return _mean_from_sum_and_count(raster)


def merge_variance(rasters: Sequence[str], mean: Raster) -> Raster:
    raster = merge(
        rasters,
        method=_copy_variance_sum_and_count(mean.pixels),
        output_count=2,
    )
    return _mean_from_sum_and_count(raster)


def merge_stddev(rasters: Sequence[str], mean: Raster) -> Raster:
    variance_raster = merge_variance(rasters, mean)
    return Raster(
        numpy.sqrt(variance_raster.pixels),
        variance_raster.transform,
        variance_raster.crs,
    )


def _mean_from_sum_and_count(raster: Raster) -> Raster:
    pixels, transform, crs = raster
    pixels_sum, pixels_count = pixels
    pixels_mean = pixels_sum / pixels_count
    return Raster(pixels_mean, transform, crs)


def _copy_sum_and_count(merged_data, new_data, merged_mask, new_mask, **kwargs):
    """
    Combines rasterio's builtin `copy_sum` and `copy_count` functions.

    Expects a 3D array of length 2, which you can get by passing
    `output_count=2` to `rasterio.merge.merge`. We split this into two arrays,
    using the first for the sum and the second for the count.
    """
    assert merged_data.ndim == 3 and len(merged_data) == 2
    merged_sum, merged_count = numpy.split(merged_data, 2)
    merged_sum_mask, merged_count_mask = numpy.split(merged_mask, 2)

    copy_sum(merged_sum, new_data, merged_sum_mask, new_mask, **kwargs)
    copy_count(merged_count, new_data, merged_count_mask, new_mask, **kwargs)


def _copy_variance_sum_and_count(mean):
    def _copy(merged_data, new_data, merged_mask, new_mask, **kwargs):
        assert merged_data.ndim == 3 and len(merged_data) == 2
        merged_sum, merged_count = numpy.split(merged_data, 2)
        merged_sum_mask, merged_count_mask = numpy.split(merged_mask, 2)

        variance = (new_data - mean) ** 2
        copy_sum(merged_sum, variance, merged_sum_mask, new_mask, **kwargs)
        copy_count(merged_count, new_data, merged_count_mask, new_mask, **kwargs)

    return _copy
