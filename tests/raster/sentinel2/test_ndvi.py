import os
from functools import wraps
from glob import glob
from tempfile import TemporaryDirectory

import boto3
import geopandas
import pytest
import rasterio.transform
import responses
from moto import mock_aws

from demeter.raster.sentinel2.constants import S3_BUCKET_NAME
from demeter.raster.sentinel2.ndvi import (
    _SAVE_TEST_FIXTURES,
    fetch_and_build_ndvi_rasters,
)

SEARCH_RESPONSES_FIXTURE_DIR = "tests/raster/fixtures/sentinel2/search_responses/"


def record_or_replay_search_responses(fixture_name):
    """
    If `_SAVE_TEST_FIXTURES` is True, records responses from the Sentinel-2
    OData search API.
    """
    if _SAVE_TEST_FIXTURES:
        from responses import _recorder

        def decorator(test_fn):
            fixture_path = os.path.join(
                SEARCH_RESPONSES_FIXTURE_DIR, f"{fixture_name}.yaml"
            )
            return _recorder.record(file_path=fixture_path)(test_fn)

        return decorator

    return replay_search_responses(fixture_name)


def replay_search_responses(fixture_name):
    """
    Replay Sentinel-2 OData API search responses, to avoid hitting the real API
    in tests.
    """

    def decorator(test_fn):
        fixture_path = os.path.join(
            SEARCH_RESPONSES_FIXTURE_DIR, f"{fixture_name}.yaml"
        )

        def wrapper(fn):
            @wraps(fn)
            @responses.activate
            def wrapped(*args, **kwargs):
                responses._add_from_file(file_path=fixture_path)
                return test_fn(*args, **kwargs)

            return wrapped

        return wrapper(test_fn)

    return decorator


@pytest.fixture(autouse=True)
def cache_directory(monkeypatch):
    if _SAVE_TEST_FIXTURES:
        yield ".sentinel2_cache"
    else:
        with TemporaryDirectory() as tmpdir:
            monkeypatch.setenv("SENTINEL2_CACHED_RASTER_FILES_DIRECTORY", tmpdir)
            yield tmpdir


@pytest.fixture(autouse=True)
def copernicus_s3_credentials(monkeypatch):
    if _SAVE_TEST_FIXTURES:
        if not os.environ.get("COPERNICUS_AWS_ENDPOINT_URL"):
            raise Exception(
                "Set COPERNICUS_AWS_ENDPOINT_URL and credentials to save test fixtures"
            )
    else:
        monkeypatch.delenv("COPERNICUS_AWS_ENDPOINT_URL", raising=False)
        monkeypatch.setenv("COPERNICUS_AWS_ACCESS_KEY_ID", "key")
        monkeypatch.setenv("COPERNICUS_AWS_SECRET_ACCESS_KEY", "secret")


@pytest.fixture(scope="module")
def copernicus_s3():
    if _SAVE_TEST_FIXTURES:
        yield None
    else:
        with mock_aws():
            s3 = boto3.client("s3")
            s3.create_bucket(Bucket=S3_BUCKET_NAME)
            yield s3


@pytest.fixture
def geometries():
    """
    Chosen to span multiple tiles.
    """
    return geopandas.read_file(
        "tests/raster/fixtures/fields_spanning_sentinel2_tiles.geojson"
    )


@pytest.fixture
def geometries_spanning_datatake_edge():
    """
    For testing detector footprint masking.
    """
    return geopandas.read_file("tests/raster/fixtures/texas_west.geojson")


@pytest.fixture(scope="module")
def sentinel2_rasters_in_s3(copernicus_s3):
    """
    These are real rasters, cropped to cover only the input geometries (plus a
    small buffer) to keep file size small.

    To regenerate these fixtures, set `_SAVE_TEST_FIXTURES` to True and rerun
    the tests.
    """
    if _SAVE_TEST_FIXTURES:
        return

    fixtures_dir = "tests/raster/fixtures/sentinel2/eodata/"

    for path in glob(os.path.join(fixtures_dir, "**/*"), recursive=True):
        _, ext = os.path.splitext(path)
        if ext not in {".safe", ".jp2"}:
            continue

        key = os.path.relpath(path, fixtures_dir)
        copernicus_s3.upload_file(path, S3_BUCKET_NAME, key)


@record_or_replay_search_responses("2024_09_fields_spanning_sentinel2_tiles")
def test_fetch_and_build_ndvi_rasters(
    geometries,
    sentinel2_rasters_in_s3,
):
    rasters = list(fetch_and_build_ndvi_rasters(geometries, 2024, 9))
    assert len(rasters) == 1  # geometries are all in UTM zone 14

    raster = rasters[0]
    assert raster.crs == "EPSG:32614"
    assert raster.mean

    assert raster.mean.shape == (1521, 319)
    assert raster.mean.pixels.count() == 12287
    assert round(raster.mean.pixels.mean(), 3) == 0.548
    assert raster.mean.crs == "EPSG:32614"

    # Check that raster bounds are within 10m (1 pixel) of input geometry bounds:
    height, width = raster.mean.shape
    raster_bounds = rasterio.transform.array_bounds(
        height, width, raster.mean.transform
    )
    input_geometry_bounds = geometries.to_crs(raster.crs).total_bounds
    assert all(abs(input_geometry_bounds - raster_bounds) < 10)


@replay_search_responses("2024_09_fields_spanning_sentinel2_tiles")
def test_fetch_and_build_ndvi_rasters_min(
    geometries,
    sentinel2_rasters_in_s3,
):
    rasters = list(
        fetch_and_build_ndvi_rasters(geometries, 2024, 9, statistics=["min"])
    )
    assert len(rasters) == 1

    raster = rasters[0]
    assert raster.min

    pixels, transform, crs = raster.min
    assert round(pixels.mean(), 3) == 0.458


@replay_search_responses("2024_09_fields_spanning_sentinel2_tiles")
def test_fetch_and_build_ndvi_rasters_max(
    geometries,
    sentinel2_rasters_in_s3,
):
    rasters = list(
        fetch_and_build_ndvi_rasters(geometries, 2024, 9, statistics=["max"])
    )
    assert len(rasters) == 1

    raster = rasters[0]
    assert raster.max

    pixels, transform, crs = raster.max
    assert round(pixels.mean(), 3) == 0.602


@replay_search_responses("2024_09_fields_spanning_sentinel2_tiles")
def test_fetch_and_build_ndvi_rasters_stddev(
    geometries,
    sentinel2_rasters_in_s3,
):
    with pytest.raises(ValueError):
        list(fetch_and_build_ndvi_rasters(geometries, 2024, 9, statistics=["stddev"]))

    rasters = list(
        fetch_and_build_ndvi_rasters(geometries, 2024, 9, statistics=["mean", "stddev"])
    )
    assert len(rasters) == 1

    raster = rasters[0]
    assert raster.mean
    assert raster.stddev

    pixels, transform, crs = raster.stddev
    assert round(pixels.mean(), 3) == 0.049


@record_or_replay_search_responses("2024_09_agoro_shea_flanagan_west")
def test_detector_footprint_mask(
    geometries_spanning_datatake_edge,
    sentinel2_rasters_in_s3,
):
    """
    Applying the detector footprint mask should prevent artifacts at the edges
    of satellite's field of view. These artifacts are most visible in the min
    NDVI raster.
    """
    rasters = list(
        fetch_and_build_ndvi_rasters(
            geometries_spanning_datatake_edge, 2024, 9, statistics=["min"]
        )
    )
    assert len(rasters) == 1

    raster = rasters[0]
    assert raster.min

    pixels, transform, crs = raster.min

    # Without applying the detector footprint mask, this raster has a min value
    # of -0.715
    assert round(pixels.min(), 3) == 0.017
