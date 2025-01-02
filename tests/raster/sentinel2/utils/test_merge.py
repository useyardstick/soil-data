import numpy
import pytest
import rasterio

from soil_data.raster.sentinel2.utils.merge import (
    merge_max,
    merge_mean,
    merge_min,
    merge_stddev,
    merge_variance,
)


@pytest.fixture
def raster_paths(tmp_path):
    raster_arrays = [
        numpy.array(
            [
                [4.0, 3.0],
                [5.0, 5.0],
            ]
        ),
        numpy.array(
            [
                [6.0, numpy.nan],
                [9.0, 4.0],
            ]
        ),
    ]

    for index, raster in enumerate(raster_arrays):
        with rasterio.open(
            tmp_path / f"raster_{index}.tif",
            "w",
            count=1,
            height=raster.shape[0],
            width=raster.shape[1],
            dtype=raster.dtype,
            nodata=numpy.nan,
            transform=rasterio.Affine(1, 0, 0, 0, -1, 0),
            crs="EPSG:4326",
        ) as dst:
            dst.write(raster, indexes=1)

    return list(tmp_path.glob("*.tif"))


def test_merge_min(raster_paths):
    min_raster, _, _ = merge_min(raster_paths)
    assert (
        min_raster
        == numpy.array(
            [
                [4.0, 3.0],
                [5.0, 4.0],
            ]
        )
    ).all()


def test_merge_max(raster_paths):
    max_raster, _, _ = merge_max(raster_paths)
    assert (
        max_raster
        == numpy.array(
            [
                [6.0, 3.0],
                [9.0, 5.0],
            ]
        )
    ).all()


def test_merge_mean(raster_paths):
    mean_raster, _, _ = merge_mean(raster_paths)
    assert (
        mean_raster
        == numpy.array(
            [
                [5.0, 3.0],
                [7.0, 4.5],
            ]
        )
    ).all()


def test_merge_variance(raster_paths):
    mean = merge_mean(raster_paths)
    variance_raster, _, _ = merge_variance(raster_paths, mean)
    assert (
        variance_raster
        == numpy.array(
            [
                [1.0, 0.0],
                [4.0, 0.25],
            ]
        )
    ).all()


def test_merge_stddev(raster_paths):
    mean = merge_mean(raster_paths)
    stddev_raster, _, _ = merge_stddev(raster_paths, mean)
    assert (
        stddev_raster
        == numpy.array(
            [
                [1.0, 0.0],
                [2.0, 0.5],
            ]
        )
    ).all()
