import os
import warnings
from collections.abc import Callable, Sequence
from contextlib import ExitStack, nullcontext
from typing import Literal, Union

import numpy
import rasterio
import rasterio.merge
from rasterio.merge import copy_count, copy_sum

from demeter.raster import Raster

MergeMethod = Union[
    Literal["first", "last", "min", "max", "sum", "count", "mean"], Callable
]


def merge(rasters: Sequence, *, method: MergeMethod = "first", **kwargs) -> Raster:
    """
    Wraps `rasterio.merge.merge` to operate on Raster instances as well as
    rasterio datasets.

    The `method` argument specifies how to handle overlapping pixels. See
    https://rasterio.readthedocs.io/en/stable/api/rasterio.merge.html for
    details on the available methods.

    In addition to rasterio's built-in methods listed above, this also supports
    a `mean` method that returns the mean of all valid overlapping pixels.
    """
    if isinstance(rasters[0], Raster):
        with ExitStack() as stack:
            datasets = [stack.enter_context(raster.as_dataset()) for raster in rasters]
            return _merge(datasets, method=method, **kwargs)

    return _merge(rasters, method=method, **kwargs)


def merge_min(rasters: Sequence, **kwargs) -> Raster:
    """
    Merge the given rasters, using the minimum value at each overlapping pixel.
    """
    return merge(rasters, method="min", **kwargs)


def merge_max(rasters: Sequence, **kwargs) -> Raster:
    """
    Merge the given rasters, using the maxiumum value at each overlapping pixel.
    """
    return merge(rasters, method="max", **kwargs)


def merge_mean(rasters: Sequence, **kwargs) -> Raster:
    """
    Merge the given rasters, using the mean value at each overlapping pixel.
    """
    return merge(rasters, method="mean", **kwargs)


def merge_variance(rasters: Sequence, mean: Raster, **kwargs) -> Raster:
    """
    Calculate the mean variance of rasters from the given mean.
    """
    raster = merge(
        rasters,
        method=_copy_variance_sum_and_count(mean.pixels),
        output_count=2,
        **kwargs,
    )
    return _mean_from_sum_and_count(raster)


def merge_stddev(rasters: Sequence, mean: Raster, **kwargs) -> Raster:
    """
    Calculate the mean standard deviation of rasters from the given mean.
    """
    variance_raster = merge_variance(rasters, mean, **kwargs)
    return Raster(
        numpy.sqrt(variance_raster.pixels),
        variance_raster.transform,
        variance_raster.crs,
    )


def _merge(sources: Sequence, *, method, output_count=None, **kwargs) -> Raster:
    # Get the CRS from the first raster. If any of the other rasters have a
    # different CRS, the call to `rasterio.merge.merge` below will raise an
    # exception, so we can safely assume this is the CRS to use for the output
    # raster.
    first_source = sources[0]

    if isinstance(first_source, (str, os.PathLike)):
        dataset_opener = rasterio.open
    else:
        dataset_opener = nullcontext

    with dataset_opener(first_source) as dataset:
        crs = dataset.crs
        num_bands = dataset.count

    if crs is None:
        raise ValueError("Rasters have no CRS")

    calculating_mean = method == "mean"
    if calculating_mean:
        if num_bands > 1:
            raise ValueError(
                "Calculating mean for multi-band rasters not yet supported"
            )

        # Stack two numpy arrays, the first with the sum of valid values at
        # each pixel, and the second with the count of valid values at each
        # pixel. Then divide the first array by the second to get the mean.
        method = _copy_sum_and_count
        output_count = 2

    pixels, transform = rasterio.merge.merge(
        sources,
        masked=True,
        method=method,
        output_count=output_count,
        **kwargs,
    )
    raster = Raster(pixels, transform, str(crs))

    if calculating_mean:
        return _mean_from_sum_and_count(raster)

    return raster


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


class OverlappingPixelsWarning(Warning):
    pass


def check_for_overlapping_pixels(
    merged_data, new_data, merged_mask, new_mask, **kwargs
):
    """
    When passed as the `method` argument to `rasterio.merge.merge`, this
    function checks whether any two rasters have data for the same pixel.
    If they do, it logs a warning.
    """
    # `merged_mask` and `new_mask` are boolean arrays with True values for
    # invalid pixels. For every pixel, one or both of `merged_mask` and
    # `new_mask` should be True. If *both* are False, it means two rasters have
    # valid data for the same pixel. If we see this, and the values are
    # different, log a warning.
    overlap_mask = ~(merged_mask | new_mask)
    if (merged_data[overlap_mask] != new_data[overlap_mask]).any():
        warnings.warn(
            "Input rasters have overlapping pixels with different values!",
            category=OverlappingPixelsWarning,
        )

    # Carry on with rasterio's default merge behavior:
    rasterio.merge.copy_first(merged_data, new_data, merged_mask, new_mask, **kwargs)
