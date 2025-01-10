from dataclasses import dataclass

import numpy
import rasterio
import rasterio.transform


@dataclass
class Raster:
    """
    Rasterio has a file-centric API. It has the concept of a `MemoryFile` for
    in-memory processing, but it's a bit clunky. This is intended as a simpler
    in-memory representation of raster data, with direct access to the raster
    pixels as a numpy masked array.

    Tools in this library return `Raster` instances. To save to disk, use the
    `save` method:

    ```python
    raster.save("path/to/file.tif")
    ```
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
            width=width,
            height=height,
            count=self.count,
            crs=self.crs,
            dtype=self.dtype,
            nodata=self.pixels.fill_value,
            transform=self.transform,
            **kwargs,
        ) as dataset:
            dataset.write(self.pixels)

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
