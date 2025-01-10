import numpy
import pytest
import rasterio

from demeter.raster import Raster


@pytest.fixture
def single_band_2d_array():
    return numpy.ma.array(
        [
            [6, 0],
            [9, 4],
        ]
    )


@pytest.fixture
def single_band_3d_array(single_band_2d_array):
    return numpy.expand_dims(single_band_2d_array, axis=0)


@pytest.fixture
def multiband_array():
    return numpy.ma.array(
        [
            [
                [6, 0],
                [9, 4],
            ],
            [
                [4, 3],
                [5, 5],
            ],
        ]
    )


@pytest.fixture
def transform():
    return rasterio.Affine(1, 0, 0, 0, -1, 0)


@pytest.fixture
def crs():
    return "EPSG:5070"


def test_load_and_save_raster(tmp_path):
    tiff_save_path = tmp_path / "testing_save.tif"

    raster = Raster.from_file("tests/raster/fixtures/polaris/lat4142_lon-88-87.tif")
    raster.save(tiff_save_path)

    raster2 = Raster.from_file(tiff_save_path)
    assert raster.transform == raster2.transform
    assert (raster.pixels == raster2.pixels).all()
    assert raster.crs == raster2.crs


def test_2d_raster(single_band_2d_array, transform, crs):
    raster = Raster(single_band_2d_array, transform, crs)
    assert raster.shape == (2, 2)
    assert raster.pixels.shape == (1, 2, 2)
    assert raster.count == 1
    assert raster.value_at(0, 1) == 9


def test_3d_raster(single_band_3d_array, transform, crs):
    raster = Raster(single_band_3d_array, transform, crs)
    assert raster.shape == (2, 2)
    assert raster.pixels.shape == (1, 2, 2)
    assert raster.count == 1
    assert raster.value_at(0, 1) == 9


def test_multiband_raster(multiband_array, transform, crs):
    raster = Raster(multiband_array, transform, crs)
    assert raster.shape == (2, 2)
    assert raster.pixels.shape == (2, 2, 2)
    assert raster.count == 2

    with pytest.raises(ValueError):
        raster.value_at(0, 1)

    assert raster.values_at(0, 1) == [9, 5]
