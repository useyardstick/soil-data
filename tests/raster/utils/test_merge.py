import numpy
import pytest
import rasterio

from demeter.raster.utils.merge import (
    merge,
    merge_max,
    merge_mean,
    merge_min,
    merge_stddev,
    merge_variance,
)


@pytest.fixture
def int_rasters_with_zero_nodata(tmp_path):
    return _save_rasters(
        tmp_path,
        [
            numpy.array(
                [
                    [6, 0],
                    [9, 4],
                ]
            ),
            numpy.array(
                [
                    [4, 3],
                    [5, 5],
                ]
            ),
        ],
        nodata=0,
    )


@pytest.fixture
def int_rasters_with_nonzero_nodata(tmp_path):
    return _save_rasters(
        tmp_path,
        [
            numpy.array(
                [
                    [6, -9999],
                    [9, 4],
                ]
            ),
            numpy.array(
                [
                    [4, 3],
                    [5, 5],
                ]
            ),
        ],
        nodata=-9999,
    )


@pytest.fixture
def float_rasters(tmp_path):
    return _save_rasters(
        tmp_path,
        [
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
        ],
        nodata=numpy.nan,
    )


def test_merge_int_rasters_with_nonzero_nodata(int_rasters_with_nonzero_nodata):
    """
    There was a bug in `rasterio.merge` that caused it to return invalid data
    when merging int rasters with a nonzero nodata value. The bug was fixed in
    rasterio 1.4.3. This is a regression test to make sure it stays fixed.
    """
    merged = merge(int_rasters_with_nonzero_nodata)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.pixels.fill_value == -9999


def test_merge_int_rasters_with_nonzero_nodata_as_float(
    int_rasters_with_nonzero_nodata,
):
    """
    The rasterio bug above also happened when converting the merged output to a
    float dtype. Make sure that stays fixed too.
    """
    merged = merge(int_rasters_with_nonzero_nodata, dtype=float)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.pixels.fill_value == -9999


def test_merge_int_rasters_with_nonzero_nodata_passing_zero_nodata(
    int_rasters_with_nonzero_nodata,
):
    """
    Passing a zero nodata value to `merge` works.
    """
    merged = merge(int_rasters_with_nonzero_nodata, nodata=0)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.pixels.fill_value == 0


def test_merge_int_rasters_with_zero_nodata(int_rasters_with_zero_nodata):
    """
    Merging int rasters with a zero nodata value works as expected.
    """
    merged = merge(int_rasters_with_zero_nodata)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.pixels.fill_value == 0


def test_merge_int_rasters_with_zero_nodata_passing_nonzero_nodata(
    int_rasters_with_zero_nodata,
):
    """
    Merging int rasters with a zero nodata value works as expected, even when
    passing a nonzero nodata value to `merge`.
    """
    merged = merge(int_rasters_with_zero_nodata, nodata=-9999)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.pixels.fill_value == -9999


def test_merge_min(float_rasters):
    min_raster, _, _ = merge_min(float_rasters)
    assert (
        min_raster
        == numpy.array(
            [
                [4.0, 3.0],
                [5.0, 4.0],
            ]
        )
    ).all()


def test_merge_max(float_rasters):
    max_raster, _, _ = merge_max(float_rasters)
    assert (
        max_raster
        == numpy.array(
            [
                [6.0, 3.0],
                [9.0, 5.0],
            ]
        )
    ).all()


def test_merge_mean(float_rasters):
    mean_raster, _, _ = merge_mean(float_rasters)
    assert (
        mean_raster
        == numpy.array(
            [
                [5.0, 3.0],
                [7.0, 4.5],
            ]
        )
    ).all()


def test_merge_variance(float_rasters):
    mean = merge_mean(float_rasters)
    variance_raster, _, _ = merge_variance(float_rasters, mean)
    assert (
        variance_raster
        == numpy.array(
            [
                [1.0, 0.0],
                [4.0, 0.25],
            ]
        )
    ).all()


def test_merge_stddev(float_rasters):
    mean = merge_mean(float_rasters)
    stddev_raster, _, _ = merge_stddev(float_rasters, mean)
    assert (
        stddev_raster
        == numpy.array(
            [
                [1.0, 0.0],
                [2.0, 0.5],
            ]
        )
    ).all()


def _save_rasters(tmp_path, arrays, nodata):
    for index, array in enumerate(arrays):
        height, width = array.shape
        with rasterio.open(
            tmp_path / f"raster_{index}.tif",
            "w",
            count=1,
            height=height,
            width=width,
            dtype=array.dtype,
            transform=rasterio.Affine(1, 0, 0, 0, -1, 0),
            crs="EPSG:4326",
            nodata=nodata,
        ) as dst:
            dst.write(array, indexes=1)

    return list(tmp_path.glob("*.tif"))
