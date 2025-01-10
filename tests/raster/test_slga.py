import os
from urllib.parse import urlparse

import geopandas
import pytest
import rasterio

from demeter.raster import slga


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
                            (117.66800871261091, -33.47764501416449),
                            (117.66800871261091, -33.48305309988765),
                            (117.68625977935301, -33.48305309988765),
                            (117.68625977935301, -33.47764501416449),
                            (117.66800871261091, -33.47764501416449),
                        ]
                    ],
                },
                "properties": {},
            }
        ]
    )


@pytest.fixture
def mock_slga(monkeypatch):
    real_open = rasterio.open

    def mock_open(url):
        path = urlparse(url).path
        filename = os.path.basename(path)
        local_path = os.path.join("tests/raster/fixtures/slga", filename)
        return real_open(local_path)

    monkeypatch.setattr("rasterio.open", mock_open)


def test_estimate_carbon_stock(mock_slga, geometries):
    rasters = slga.estimate_carbon_stock(
        geometries,
        end_depth=30,
    )
    assert rasters.stddev
    assert rasters.mean.shape == rasters.stddev.shape == (8, 23)
    assert rasters.mean.transform == rasters.stddev.transform
    assert rasters.mean.crs == "EPSG:4326"
    assert round(rasters.mean.pixels.mean(), 3) == 1.843
    assert round(rasters.stddev.pixels.mean(), 3) == 1.169


def test_fetch_slga_data_for_depth_range(mock_slga, geometries):
    rasters = slga.fetch_slga_data_for_depth_range(
        geometries,
        slga.SoilProperty.PH,
        end_depth=30,
    )
    assert rasters.stddev
    assert rasters.mean.shape == rasters.stddev.shape == (8, 23)
    assert rasters.mean.transform == rasters.stddev.transform
    assert rasters.mean.crs == "EPSG:4326"
    assert round(rasters.mean.pixels.mean(), 3) == 6.212
    assert round(rasters.stddev.pixels.mean(), 3) == 0.806
