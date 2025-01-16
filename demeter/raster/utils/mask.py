import geopandas
import rasterio.mask

from demeter.raster import Raster


# TODO: combine these
def mask(dataset, shapes, **kwargs) -> Raster:
    """
    Wraps rasterio.mask.mask to return a Raster instance instead of a
    (raster, transform) 2-tuple.
    """
    crs = dataset.crs
    assert crs

    pixels, transform = rasterio.mask.mask(dataset, shapes, filled=False, **kwargs)
    return Raster(pixels, transform, str(crs))


def mask_raster(raster: Raster, shapes, **kwargs) -> Raster:
    """
    Wraps `rasterio.mask.mask` to operate on in-memory rasters instead of
    rasterio datasets.
    """
    if isinstance(shapes, geopandas.GeoDataFrame):
        shapes = shapes.geometry

    with raster.as_dataset() as dataset:
        return mask(dataset, shapes, **kwargs)
