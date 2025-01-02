import geopandas
import numpy
import pytest
import rasterio.transform

from soil_data.raster.usgs.topography import RASTER_CRS, fetch_and_merge_rasters


@pytest.fixture
def geometries():
    # Chosen to span 2 USGS elevation tiles:
    return geopandas.read_file("tests/raster/fixtures/permian_basin.geojson")


# TODO: save test fixtures
def test_fetch_and_merge_rasters(geometries):
    raster, transform, crs = fetch_and_merge_rasters(geometries)
    assert raster.shape == (14934, 4789)
    assert raster.count() == 6008
    assert crs == "EPSG:4269"

    # Check that raster bounds are within 5m of input geometry bounds:
    height, width = raster.shape
    raster_bounds = rasterio.transform.array_bounds(height, width, transform)
    input_geometry_bounds = geometries.to_crs(RASTER_CRS).total_bounds
    assert all(abs(input_geometry_bounds - raster_bounds) < 5)
