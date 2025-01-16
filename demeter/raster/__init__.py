from dataclasses import dataclass
from typing import Any

import numpy
import rasterio
import rasterio.transform

from demeter.raster.utils.transform import extract_resolution_from_transform


@dataclass
class Raster:
    """
    Rasterio has a file-centric API. It has the concept of a MemoryFile for
    in-memory processing, but it's a bit clunky. This is intended as a simpler
    in-memory representation of raster data, with direct access to the raster
    pixels as a numpy masked array.
    """

    pixels: numpy.ma.MaskedArray
    """
    A 3-dimensional array: one 2D array per band.
    """

    transform: rasterio.Affine
    """
    The transform from pixel coordinates to geographic coordinates in this
    raster's CRS.
    """

    crs: str
    """
    This raster's CRS, as a string. For example: "EPSG:5070".
    """

    @classmethod
    def from_file(cls, path: str) -> "Raster":
        """
        Read the file at the given path into memory and return a Raster
        instance.
        """
        with rasterio.open(path) as dataset:
            transform = dataset.transform
            crs = dataset.crs
            pixels = dataset.read(masked=True)

        return cls(pixels, transform, str(crs))

    @property
    def count(self) -> int:
        """
        Number of bands.
        """
        return self.pixels.shape[0]

    @property
    def shape(self) -> tuple[int, int]:
        """
        Height and width of the raster.
        """
        return self.pixels.shape[-2:]

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """
        The raster's bounds, in the raster's CRS.
        """
        height, width = self.shape
        return rasterio.transform.array_bounds(height, width, self.transform)

    @property
    def resolution(self) -> tuple[float, float]:
        """
        The raster's (x, y) resolution.
        """
        return extract_resolution_from_transform(self.transform)

    @property
    def grid_offset(self) -> tuple[float, float]:
        """
        The (x, y) offset of this raster's origin point on a grid aligned with
        its resolution.
        """
        xres, yres = self.resolution
        return self.transform.xoff % xres, self.transform.yoff % yres

    @property
    def dtype(self):
        """
        The raster's data type.
        """
        return self.pixels.dtype

    def value_at(self, x: float, y: float):
        """
        Find the pixel corresponding to the given coordinates, and return its
        value. Assumes a single-band raster.
        """
        if self.count > 1:
            raise ValueError("Raster has multiple bands. Use `values_at` instead.")

        return self.values_at(x, y)[0]

    def values_at(self, x: float, y: float) -> list:
        """
        Find all pixels at the given coordinates from all bands, and return
        them as a list.
        """
        row, col = rasterio.transform.rowcol(self.transform, x, y)
        return self.pixels[:, int(row), int(col)].tolist()

    def save(self, path: str, **kwargs):
        """
        Save the raster to the given path.
        """
        height, width = self.shape
        with rasterio.open(
            path,
            mode="w",
            **self._rasterio_profile(),
            **kwargs,
        ) as dataset:
            dataset.write(self.pixels)

    # TODO: add tempfile implementation
    def as_dataset(self) -> rasterio.io.DatasetReader:
        """
        Write this raster to a `MemoryFile`, then open it for reading.
        """
        memory_file = rasterio.io.MemoryFile()

        height, width = self.shape
        with memory_file.open(driver="GTiff", **self._rasterio_profile()) as dataset:
            dataset.write(self.pixels)

        # Open a second time for reading:
        dataset = memory_file.open()

        # To ensure the memory file is closed when the dataset is closed, add it to
        # the dataset's exit stack. This is similar to what rasterio does here:
        # https://github.com/rasterio/rasterio/blob/1.4.3/rasterio/__init__.py#L301
        dataset._env.enter_context(memory_file)
        return dataset

    def _rasterio_profile(self) -> dict[str, Any]:
        height, width = self.shape
        return {
            "width": width,
            "height": height,
            "count": self.count,
            "crs": self.crs,
            "dtype": self.dtype,
            "nodata": self.pixels.fill_value,
            "transform": self.transform,
        }

    def __post_init__(self):
        """
        Runtime validation to ensure:

        - All attributes are set, in case mypy doesn't catch a missing value.
        - The `pixels` array is a 3-dimensional masked array.
        """
        if self.pixels is None:
            raise ValueError("Raster has no pixels")
        if not self.transform:
            raise ValueError("Raster has no transform")
        if not self.crs:
            raise ValueError("Raster has no CRS")

        if not isinstance(self.pixels, numpy.ma.MaskedArray):
            raise ValueError("Raster pixels must be a masked array")
        if self.pixels.ndim == 2:
            self.pixels = numpy.expand_dims(self.pixels, axis=0)
        elif self.pixels.ndim != 3:
            raise ValueError("Raster array must be 2 or 3-dimensional")

    def __iter__(self):
        """
        Raster used to be a NamedTuple. This provides backward-compatibility
        for unpacking, as in:

        ```
        pixels, transform, crs = raster
        ```
        """
        return iter((self.pixels, self.transform, self.crs))
