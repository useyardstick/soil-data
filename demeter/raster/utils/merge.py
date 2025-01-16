import os
import warnings
from collections.abc import Sequence
from contextlib import ExitStack, nullcontext

import rasterio
import rasterio.merge

from demeter.raster import Raster


def merge(rasters: Sequence, **kwargs) -> Raster:
    """
    Wraps `rasterio.merge.merge` to operate on Raster instances as well as
    rasterio datasets.
    """
    if isinstance(rasters[0], Raster):
        with ExitStack() as stack:
            datasets = [stack.enter_context(raster.as_dataset()) for raster in rasters]
            return _merge(datasets, **kwargs)

    return _merge(rasters, **kwargs)


def _merge(sources: Sequence, **kwargs) -> Raster:
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

    if crs is None:
        raise ValueError("Rasters have no CRS")

    pixels, transform = rasterio.merge.merge(sources, masked=True, **kwargs)
    return Raster(pixels, transform, str(crs))


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
