import geopandas
import numpy
import pytest
import rasterio

from demeter.raster.utils import Raster, mask_raster, merge


def test_mask_raster():
    matrix = numpy.ma.ones((4, 4))
    shapes = geopandas.GeoDataFrame.from_features(
        [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            (1, 1),
                            (3, 1),
                            (3, 3),
                            (1, 3),
                            (1, 1),
                        ]
                    ],
                },
                "properties": {},
            },
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            (2, 2),
                            (4, 2),
                            (4, 4),
                            (2, 4),
                            (2, 2),
                        ]
                    ],
                },
                "properties": {},
            },
        ]
    )
    result, *_ = mask_raster(
        Raster(matrix, transform=rasterio.Affine.identity(), crs="EPSG:4326"),
        shapes=shapes,
    )
    expected = numpy.ma.array(
        matrix,
        mask=~numpy.ma.make_mask(
            [
                [0, 0, 0, 0],
                [0, 1, 1, 0],
                [0, 1, 1, 1],
                [0, 0, 1, 1],
            ]
        ),
    )
    assert numpy.ma.allequal(result, expected)


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


def test_merge_int_rasters_with_nonzero_nodata(int_rasters_with_nonzero_nodata):
    """
    There was a bug in `rasterio.merge` that caused it to return invalid data
    when merging int rasters with a nonzero nodata value. The bug was fixed in
    rasterio 1.4.3. This is a regression test to make sure it stays fixed.
    """
    merged = merge(int_rasters_with_nonzero_nodata)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.nodata == -9999


def test_merge_int_rasters_with_nonzero_nodata_as_float(
    int_rasters_with_nonzero_nodata,
):
    """
    The rasterio bug above also happened when converting the merged output to a
    float dtype. Make sure that stays fixed too.
    """
    merged = merge(int_rasters_with_nonzero_nodata, dtype=float)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.nodata == -9999


def test_merge_int_rasters_with_nonzero_nodata_passing_zero_nodata(
    int_rasters_with_nonzero_nodata,
):
    """
    Passing a zero nodata value to `merge` works.
    """
    merged = merge(int_rasters_with_nonzero_nodata, nodata=0)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.nodata == 0


def test_merge_int_rasters_with_zero_nodata(int_rasters_with_zero_nodata):
    """
    Merging int rasters with a zero nodata value works as expected.
    """
    merged = merge(int_rasters_with_zero_nodata)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.nodata == 0


def test_merge_int_rasters_with_zero_nodata_passing_nonzero_nodata(
    int_rasters_with_zero_nodata,
):
    """
    Merging int rasters with a zero nodata value works as expected, even when
    passing a nonzero nodata value to `merge`.
    """
    merged = merge(int_rasters_with_zero_nodata, nodata=-9999)
    assert numpy.ma.allequal(merged.pixels, numpy.ma.masked_array([[6, 3], [9, 4]]))
    assert merged.nodata == -9999


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
