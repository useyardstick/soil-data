from contextlib import nullcontext

import geopandas
import rasterio.mask

from demeter.raster import Raster


def mask(raster, shapes, **kwargs) -> Raster:
    """
    Wraps `rasterio.mask.mask`, with the following differences:

    - Can accept a `Raster` instance as well as a rasterio dataset.
    - Returns a `Raster` instance instead of a (raster, transform) 2-tuple.
    """
    if isinstance(shapes, geopandas.GeoDataFrame):
        shapes = shapes.geometry

    if isinstance(raster, Raster):
        dataset_opener = raster.as_dataset
    else:
        dataset_opener = lambda: nullcontext(raster)  # type: ignore

    with dataset_opener() as dataset:
        crs = dataset.crs
        if crs is None:
            raise ValueError("Raster has no CRS")

        pixels, transform = rasterio.mask.mask(dataset, shapes, filled=False, **kwargs)

    return Raster(pixels, transform, str(crs))
