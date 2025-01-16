import math
from collections.abc import Callable, Iterable
from contextlib import ExitStack
from tempfile import NamedTemporaryFile
from typing import Literal, Optional, Union

import numpy
import rasterio
import rasterio.warp

from demeter.raster import Raster, extract_resolution_from_transform
from demeter.raster.utils.merge import merge


def reproject(
    raster: Union[str, Raster],
    crs: str,
    resampling_method: Literal[
        "nearest",
        "bilinear",
        "cubic",
        "cubic_spline",
        "lanczos",
        "average",
        "mode",
        "gauss",
        "max",
        "min",
        "med",
        "q1",
        "q3",
        "sum",
        "rms",
    ],
    align_to_transform: Optional[rasterio.Affine] = None,
) -> Raster:
    """
    Reproject a raster to a different coordinate reference system.

    This is a lossy operation as it involves resampling. The resampling
    algorithm to use depends on the nature of the data: `nearest` is a good
    choice for categorical data, whereas `bilinear` or `average` might make
    more sense for continuous data. See
    https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling
    for the list of available resampling algorithms.

    If `align_to_transform` is provided, align the output raster to the given
    transform's pixel grid. This will perform up/downsampling if the given
    transform's resolution doesn't match the input raster's resolution.
    """
    if isinstance(raster, str):
        raster = Raster.from_file(raster)

    if (
        raster.crs == crs
        and align_to_transform is None
        or align_to_transform == raster.transform
    ):
        return raster

    if align_to_transform is None:
        destination = None
        dst_transform = None
    else:
        # Calculate output transform and array shape:
        src_height, src_width = raster.shape
        (
            dst_transform,
            dst_width,
            dst_height,
        ) = rasterio.warp.calculate_default_transform(
            raster.crs,
            crs,
            src_width,
            src_height,
            *raster.bounds,
            resolution=extract_resolution_from_transform(align_to_transform),
        )

        # Initialize output array. We need this because
        # `rasterio.warp.reproject` requires a destination array when passing a
        # transform.
        destination = numpy.empty(
            (raster.count, dst_height, dst_width), dtype=raster.dtype
        )

        # Shift the transform to align with the given transform's pixel grid:
        dst_transform = _align_transform(dst_transform, align_to_transform)

    pixels, transform = rasterio.warp.reproject(
        raster.pixels,
        destination,
        src_crs=raster.crs,
        src_transform=raster.transform,
        src_nodata=raster.pixels.fill_value,
        dst_crs=crs,
        dst_transform=dst_transform,
        resampling=rasterio.enums.Resampling[resampling_method],
        masked=True,
    )

    # `rasterio.warp.reproject` doesn't return a masked array, even with
    # `masked=True`. See https://github.com/rasterio/rasterio/pull/3289
    pixels = numpy.ma.masked_equal(pixels, raster.pixels.fill_value)

    return Raster(pixels, transform, crs)


def align(
    raster: Union[str, Raster],
    to: Union[str, Raster],
    resampling_method: Literal[
        "nearest",
        "bilinear",
        "cubic",
        "cubic_spline",
        "lanczos",
        "average",
        "mode",
        "gauss",
        "max",
        "min",
        "med",
        "q1",
        "q3",
        "sum",
        "rms",
    ],
) -> Raster:
    """
    Align a raster to another raster's grid.
    """
    transform, crs = _extract_transform_and_crs(to)
    return reproject(raster, crs, resampling_method, transform)


def reproject_and_merge(
    rasters: Iterable[Union[str, Raster]],
    crs: str,
    resampling_method: Literal[
        "nearest",
        "bilinear",
        "cubic",
        "cubic_spline",
        "lanczos",
        "average",
        "mode",
        "gauss",
        "max",
        "min",
        "med",
        "q1",
        "q3",
        "sum",
        "rms",
    ],
    merge_method: Union[
        Literal["first", "last", "min", "max", "sum", "count", "mean"], Callable
    ] = "first",
    align_to_transform: Optional[rasterio.Affine] = None,
    **kwargs,
) -> Raster:
    """
    Reproject multiple rasters to a common CRS, then merge them.

    The `merge_method` argument specifies how to handle overlapping
    pixels. Other arguments are passed to `merge`. Example:

    ```python
    merged = reproject_and_merge(
        rasters,
        crs="EPSG:4326",
        resampling_method="average",  # how to resample when reprojecting
        merge_method="mean",          # how to merge overlapping pixels
    )
    ```
    """
    with ExitStack() as stack:
        paths_to_merge = []
        for raster in rasters:
            tempfile = stack.enter_context(NamedTemporaryFile(suffix=".tif"))
            reprojected_source = reproject(
                raster, crs, resampling_method, align_to_transform
            )
            reprojected_source.save(tempfile.name)
            paths_to_merge.append(tempfile.name)

        if not paths_to_merge:
            raise ValueError("No rasters to merge")

        return merge(
            paths_to_merge,
            resampling=rasterio.enums.Resampling[resampling_method],
            method=merge_method,
            **kwargs,
        )


def align_and_merge(
    rasters: Iterable[Union[str, Raster]],
    to: Union[str, Raster],
    resampling_method: Literal[
        "nearest",
        "bilinear",
        "cubic",
        "cubic_spline",
        "lanczos",
        "average",
        "mode",
        "gauss",
        "max",
        "min",
        "med",
        "q1",
        "q3",
        "sum",
        "rms",
    ],
    merge_method: Union[
        Literal["first", "last", "min", "max", "sum", "count", "mean"], Callable
    ] = "first",
    **kwargs,
) -> Raster:
    """
    Align multiple rasters to the given raster's grid, then merge them.

    Keyword arguments are passed to `merge`.
    """
    transform, crs = _extract_transform_and_crs(to)
    return reproject_and_merge(
        rasters, crs, resampling_method, merge_method, transform, **kwargs
    )


def _extract_transform_and_crs(
    raster: Union[str, Raster]
) -> tuple[rasterio.Affine, str]:
    if isinstance(raster, str):
        with rasterio.open(raster) as dataset:
            return dataset.transform, str(dataset.crs)

    return raster.transform, raster.crs


def _align_transform(
    transform: rasterio.Affine, to: rasterio.Affine
) -> rasterio.Affine:
    resolution = extract_resolution_from_transform(transform)
    if resolution != extract_resolution_from_transform(to):
        raise ValueError("Transforms must have the same resolution")

    xdiff = to.xoff - transform.xoff
    ydiff = to.yoff - transform.yoff
    xres, yres = resolution
    xoff = _calculate_min_offset(xdiff, xres)
    yoff = _calculate_min_offset(ydiff, yres)
    a, b, c, d, e, f, *rest = transform
    return rasterio.Affine(a, b, c + xoff, d, e, f + yoff, *rest)


def _calculate_min_offset(distance: float, resolution: float) -> float:
    offset = distance % math.copysign(resolution, distance)
    if abs(offset) > resolution / 2:
        offset -= math.copysign(resolution, offset)
    assert abs(offset) <= resolution / 2
    return offset
