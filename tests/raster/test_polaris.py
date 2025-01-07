import re
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import geopandas
import pytest
import rasterio.transform

from demeter.raster import polaris

FIXTURES = {
    filename: Path("tests", "raster", "fixtures", "polaris", filename).read_bytes()
    for filename in (
        "lat4142_lon-88-87.tif",
        "lat4243_lon-88-87.tif",
    )
}


@pytest.fixture(autouse=True)
def local_cache(monkeypatch):
    with TemporaryDirectory() as tmpdir:
        monkeypatch.setenv("POLARIS_CACHED_RASTER_FILES_DIRECTORY", tmpdir)
        yield tmpdir


@pytest.fixture
def remote_cache(monkeypatch, s3):
    bucket_name = "polaris-cache"
    s3.create_bucket(Bucket=bucket_name)
    monkeypatch.setenv("POLARIS_REMOTE_CACHE", f"s3://{bucket_name}")


@pytest.fixture
def geometries():
    return geopandas.GeoDataFrame.from_features(
        [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            (-88, 43),
                            (-88, 41.5),
                            (-87, 41.5),
                            (-87, 43),
                            (-88, 43),
                        ]
                    ],
                },
                "properties": {},
            }
        ]
    )


@pytest.fixture
def mock_polaris(requests_mock):
    soil_properties_pattern = "|".join(
        soil_property.value for soil_property in polaris.SoilProperty
    )
    statistics_pattern = "|".join(stat.value for stat in polaris.Statistic)
    depth_values = [depth.value for depth in polaris.Depth]
    depths_pattern = "|".join(
        f"{start_depth}_{end_depth}" for start_depth, end_depth in depth_values
    )
    return [
        requests_mock.get(
            re.compile(
                f"http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/({soil_properties_pattern})/({statistics_pattern})/({depths_pattern})/{filename}"
            ),
            content=content,
            headers={"Content-Type": "image/tiff"},
        )
        for filename, content in FIXTURES.items()
    ]


def test_fetch_polaris_data_for_depth_range(mock_polaris, geometries):
    rasters = polaris.fetch_polaris_data_for_depth_range(
        geometries,
        soil_property="om",
        end_depth=100,
    )
    assert rasters.stddev
    assert rasters.mean.pixels.shape == rasters.stddev.pixels.shape
    assert rasters.mean.transform == rasters.stddev.transform
    assert rasters.mean.crs == "EPSG:4326"


def test_fetch_polaris_data_for_depth_range_below_ground(mock_polaris, geometries):
    rasters = polaris.fetch_polaris_data_for_depth_range(
        geometries,
        soil_property="om",
        start_depth=30,
        end_depth=100,
    )
    assert rasters.stddev
    assert rasters.mean.pixels.shape == rasters.stddev.pixels.shape
    assert rasters.mean.transform == rasters.stddev.transform

    # Check that we didn't request any POLARIS tiles above 30cm:
    for fixture in mock_polaris:
        for request in fixture.request_history:
            for depth in polaris.Depth.select_between(0, 30):
                start_depth, end_depth = depth.value
                assert f"{start_depth}_{end_depth}" not in request.path


def test_fetch_polaris_data_for_depth_range_arbitrary_depths(mock_polaris, geometries):
    rasters = polaris.fetch_polaris_data_for_depth_range(
        geometries,
        soil_property="om",
        start_depth=10,
        end_depth=45,
    )
    assert rasters.stddev
    assert rasters.mean.pixels.shape == rasters.stddev.pixels.shape
    assert rasters.mean.transform == rasters.stddev.transform

    # Check that we only fetched the necessary POLARIS depths:
    requested_depths = set()
    for fixture in mock_polaris:
        for request in fixture.request_history:
            match = re.search(r"/(\d+)_(\d+)/[^/]+\.tif", request.path)
            assert match
            requested_depths.add(match.groups())

    assert requested_depths == {("5", "15"), ("15", "30"), ("30", "60")}


def test_fetch_polaris_data_for_depth_range_median_and_mode(mock_polaris, geometries):
    rasters = polaris.fetch_polaris_data_for_depth_range(
        geometries,
        soil_property=polaris.SoilProperty.ORGANIC_MATTER,
        additional_statistics=[polaris.Statistic.MEDIAN, polaris.Statistic.MODE],
        end_depth=100,
    )
    assert rasters.stddev
    assert rasters.median
    assert rasters.mode
    assert (
        rasters.mean.pixels.shape
        == rasters.stddev.pixels.shape
        == rasters.median.pixels.shape
        == rasters.mode.pixels.shape
    )
    assert (
        rasters.mean.transform
        == rasters.stddev.transform
        == rasters.median.transform
        == rasters.mode.transform
    )


def test_fetch_polaris_data(mock_polaris, geometries):
    raster, transform, *_ = polaris.fetch_polaris_data(
        geometries,
        soil_property=polaris.SoilProperty.ORGANIC_MATTER,
        statistic=polaris.Statistic.MEAN,
        depth=polaris.Depth.ZERO_TO_FIVE_CM,
    )

    # Check that the transform maps geographic coordinates to raster indices:
    height, width = raster.shape
    longitude_min, latitude_min, longitude_max, latitude_max = geometries.total_bounds
    northest, westest = rasterio.transform.rowcol(
        transform, longitude_min, latitude_max
    )
    southest, eastest = rasterio.transform.rowcol(
        transform, longitude_max, latitude_min
    )
    assert westest == 0
    assert northest == 0
    assert eastest == pytest.approx(width, abs=1)
    assert southest == pytest.approx(height, abs=1)


def test_fetch_polaris_data_with_remote_cache(
    mock_polaris, geometries, local_cache, remote_cache
):
    for _ in range(2):
        polaris.fetch_polaris_data(
            geometries,
            soil_property=polaris.SoilProperty.ORGANIC_MATTER,
            statistic=polaris.Statistic.MEAN,
            depth=polaris.Depth.ZERO_TO_FIVE_CM,
        )
        shutil.rmtree(local_cache)  # clear local cache

    assert all(request.call_count == 1 for request in mock_polaris)


def test_select_depths_for_polaris():
    depths = polaris.Depth.select_between(0, 30)
    assert {x.value for x in depths} == {(0, 5), (5, 15), (15, 30)}
    assert sum(x.thickness for x in depths) == 30

    depths = polaris.Depth.select_between(0, 100)
    assert {x.value for x in depths} == {
        (0, 5),
        (5, 15),
        (15, 30),
        (30, 60),
        (60, 100),
    }

    assert sum(x.thickness for x in depths) == 100

    depths = polaris.Depth.select_between(0, 200)
    assert {x.value for x in depths} == {
        (0, 5),
        (5, 15),
        (15, 30),
        (30, 60),
        (60, 100),
        (100, 200),
    }

    assert sum(x.thickness for x in depths) == 200

    with pytest.raises(Exception):
        polaris.Depth.select_between(0, 99)


def test_estimate_carbon_stock(mock_polaris, geometries):
    rasters = polaris.estimate_carbon_stock(
        geometries,
        end_depth=100,
    )
    assert rasters.stddev
    assert rasters.mean.pixels.shape == rasters.stddev.pixels.shape
    assert rasters.mean.transform == rasters.stddev.transform
