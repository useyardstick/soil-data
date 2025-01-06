import geopandas
import numpy
import pytest
import rasterio.transform

from demeter.raster.usgs.hydrography import (
    RASTER_CRS,
    fetch_and_merge_rasters,
    find_hu4_codes,
)


@pytest.fixture
def geometries():
    return geopandas.read_file("tests/raster/fixtures/tango_oscar.geojson")


def test_find_hu4_codes(geometries):
    assert sorted(find_hu4_codes(geometries)) == ["1022", "1023"]


# TODO: save test fixtures
def test_fetch_and_merge_rasters(geometries):
    raster, transform, crs = fetch_and_merge_rasters("cat.tif", geometries)
    assert raster.shape == (1632, 8700)
    assert raster.count() == 41196
    assert crs in {"EPSG:5070", "ESRI:102039"}
    assert numpy.array_equal(
        numpy.unique(raster).compressed(),
        [
            4487,
            10337,
            11215,
            11248,
            12413,
            12461,
            12479,
            12816,
            16420,
            18589,
            21115,
            29932,
            30755,
            34945,
        ],
    )

    # Check that raster bounds are within 5m of input geometry bounds:
    height, width = raster.shape
    raster_bounds = rasterio.transform.array_bounds(height, width, transform)
    input_geometry_bounds = geometries.to_crs(RASTER_CRS).total_bounds
    assert all(abs(input_geometry_bounds - raster_bounds) < 5)
