import geopandas
import pytest

from demeter.raster.sentinel2.utils.tiles import find_tiles_for_geometries


@pytest.fixture
def fields_spanning_tile_boundary():
    return geopandas.read_file(
        "tests/raster/fixtures/fields_spanning_sentinel2_tiles.geojson"
    )


@pytest.fixture
def fields_in_separate_tiles():
    return geopandas.read_file("tests/raster/fixtures/texas.geojson")


@pytest.fixture
def field_with_orbit_that_crosses_tile_but_not_geometry():
    return geopandas.read_file("tests/raster/fixtures/humboldt_park.geojson")


def test_find_tiles_for_geometries_on_tile_boundary(fields_spanning_tile_boundary):
    """
    One of the input geometries intersects with two tiles (14TML and 14TMM),
    but the other only intersects with 14TMM. In such cases, we want to
    download the minimum number of tiles that cover all the input geometries.
    """
    tiles = list(find_tiles_for_geometries(fields_spanning_tile_boundary))
    tile_ids = {tile.tile_id for tile in tiles}
    assert tile_ids == {"14TMM"}


def test_find_tiles_for_geometries_in_separate_tiles(fields_in_separate_tiles):
    tiles = list(find_tiles_for_geometries(fields_in_separate_tiles))
    tile_ids = {tile.tile_id for tile in tiles}
    assert tile_ids == {"13SFA", "14SKF"}


def test_find_tiles_for_geometries_spanning_tiles(fields_in_separate_tiles):
    bounding_box = fields_in_separate_tiles.dissolve().envelope
    tiles = list(find_tiles_for_geometries(bounding_box))
    tile_ids = {tile.tile_id for tile in tiles}
    assert tile_ids == {"13SFA", "13SGA", "14SKF"}


def test_find_tiles_for_geometries_excludes_orbits(
    field_with_orbit_that_crosses_tile_but_not_geometry,
):
    """
    The Sentinel-2 satellites fly over this tile on several orbits, but not
    all of those orbits cover the input geometries. Check that we only include
    the orbits that we need.
    """
    tiles = list(
        find_tiles_for_geometries(field_with_orbit_that_crosses_tile_but_not_geometry)
    )
    assert {tile.tile_id for tile in tiles} == {"16TDM"}

    # Orbit number 83 intersects with the tile, but not with the geometry.
    # Make sure it's not included:
    assert {tile.relative_orbit_number for tile in tiles} == {47, 126}
