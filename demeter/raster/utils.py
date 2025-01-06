import os
import warnings
from collections.abc import Iterable, Sequence
from contextlib import ExitStack, nullcontext
from enum import Enum
from typing import NamedTuple, TypeVar, Union

import geopandas
import numpy
import pandas
import rasterio
import rasterio.mask
import rasterio.merge
import rasterio.transform


# TODO: consider using a dataclass instead of NamedTuple so we can use
# __post_init__ for runtime type-checking
class Raster(NamedTuple):
    """
    Rasterio has a file-centric API. It has the concept of a MemoryFile for
    in-memory processing, but it's a bit clunky. This is intended as a simpler
    in-memory representation of raster data.
    """

    pixels: numpy.ma.MaskedArray
    transform: rasterio.Affine
    crs: str

    @classmethod
    def from_file(cls, path: str) -> "Raster":
        with rasterio.open(path) as dataset:
            transform = dataset.transform
            crs = dataset.crs
            assert crs
            pixels = dataset.read(1, masked=True)
            assert pixels.ndim == 2

        return cls(pixels, transform, str(crs))

    @property
    def shape(self) -> tuple[int, int]:
        return self.pixels.shape

    @property
    def dtype(self):
        return self.pixels.dtype

    @property
    def nodata(self):
        fill_value = self.pixels.fill_value
        if fill_value == numpy.ma.default_fill_value(self.dtype):
            return None

        return fill_value

    def value_at(self, x: float, y: float):
        """
        Find the pixel corresponding to the given coordinates, and return its value.
        """
        row, col = rasterio.transform.rowcol(self.transform, x, y)
        return self.pixels[int(row), int(col)]

    def save(self, path: str, masked: bool = True, **kwargs):
        height, width = self.pixels.shape
        with rasterio.open(
            path,
            mode="w",
            width=width,
            height=height,
            count=1,
            crs=self.crs,
            dtype=self.dtype,
            nodata=self.nodata,
            transform=self.transform,
            **kwargs,
        ) as dataset:
            dataset.write(self.pixels, 1, masked=masked)


def mask(dataset, shapes, **kwargs) -> Raster:
    """
    Wraps rasterio.mask.mask to return a Raster instance instead of a
    (raster, transform) 2-tuple.
    """
    assert len(dataset.shape) == 2

    crs = dataset.crs
    assert crs

    pixels, transform = rasterio.mask.mask(
        dataset, shapes, filled=False, indexes=1, **kwargs
    )
    return Raster(pixels, transform, str(crs))


def mask_raster(raster: Raster, shapes, **kwargs) -> Raster:
    """
    Wraps `rasterio.mask.mask` to operate on in-memory rasters instead of
    rasterio datasets.
    """
    if isinstance(shapes, geopandas.GeoDataFrame):
        shapes = shapes.geometry

    with raster_as_dataset(raster) as dataset:
        return mask(dataset, shapes, **kwargs)


def merge(sources: Sequence, **kwargs) -> Raster:
    """
    Wraps `rasterio.merge.merge` to return a Raster instance instead of a
    (raster, transform) 2-tuple.
    """
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
    return Raster(pixels.squeeze(), transform, str(crs))


def merge_rasters(rasters: Iterable[Raster], **kwargs) -> Raster:
    """
    Wraps `rasterio.merge.merge` to operate on in-memory rasters instead of
    rasterio datasets.
    """
    with ExitStack() as stack:
        datasets = [
            stack.enter_context(raster_as_dataset(raster)) for raster in rasters
        ]
        merged = merge(datasets, **kwargs)

    return merged


def raster_as_dataset(raster: Raster) -> rasterio.io.DatasetReader:
    memory_file = rasterio.io.MemoryFile()

    height, width = raster.shape
    with memory_file.open(
        driver="GTiff",
        count=1,
        height=height,
        width=width,
        dtype=raster.dtype,
        transform=raster.transform,
        crs=raster.crs,
        nodata=raster.nodata,
    ) as dataset:
        dataset.write(raster.pixels, 1, masked=True)

    # Open a second time for reading:
    dataset = memory_file.open()

    # To ensure the memory file is closed when the dataset is closed, add it to
    # the dataset's exit stack. This is similar to what rasterio does here:
    # https://github.com/rasterio/rasterio/blob/1.4.3/rasterio/__init__.py#L301
    dataset._env.enter_context(memory_file)
    return dataset


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


# TODO: replace this with typing.Self after upgrading to python 3.11
_DepthType = TypeVar("_DepthType", bound="DepthEnum")


class DepthEnum(Enum):
    """
    Provides helpful methods for enumerating depth ranges.

    Values should be (start_depth, end_depth) 2-tuples.
    """

    @classmethod
    def select_between(
        cls: type[_DepthType], start_depth: int, end_depth: int
    ) -> list[_DepthType]:
        start_depths = {depth.start_depth for depth in cls}
        if start_depth not in start_depths:
            raise Exception(f"start_depth {start_depth} must be one of {start_depths}")

        end_depths = {depth.end_depth for depth in cls}
        if end_depth not in end_depths:
            raise Exception(f"end_depth {end_depth} must be one of {end_depths}")

        return cls.select_including(start_depth, end_depth)

    @classmethod
    def select_including(
        cls: type[_DepthType], start_depth: int, end_depth: int
    ) -> list[_DepthType]:
        max_depth = max(depth.end_depth for depth in cls)
        if start_depth < 0 or end_depth > max_depth:
            raise Exception(f"Maximum depth range: 0 - {max_depth}")

        if end_depth <= start_depth:
            raise Exception(
                f"end_depth {end_depth} must be greater than start_depth {start_depth}"
            )

        selected_depths = (
            depth
            for depth in cls
            if start_depth < depth.end_depth and end_depth > depth.start_depth
        )
        return sorted(selected_depths, key=lambda depth: depth.start_depth)

    @property
    def start_depth(self) -> int:
        start_depth, _ = self.value
        return start_depth

    @property
    def end_depth(self) -> int:
        _, end_depth = self.value
        return end_depth

    @property
    def thickness(self) -> int:
        start_depth, end_depth = self.value
        return end_depth - start_depth


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
