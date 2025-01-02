import re
from collections.abc import Collection, Hashable, Iterable, Sequence, Set
from itertools import chain, combinations
from typing import NamedTuple, Union

import geopandas
import shapely
from pyproj.aoi import AreaOfInterest
from pyproj.database import CRSInfo, query_utm_crs_info

from soil_data.raster.sentinel2.constants import (
    ORBITS_CONUS_PATH,
    SWATH_WIDTH_KM,
    TILES_CONUS_PATH,
)
from soil_data.raster.utils import bounds_snapped_to_grid


class TileMetadata(NamedTuple):
    tile_id: str
    relative_orbit_number: int


def find_tiles_for_geometries(
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries],
) -> Iterable[TileMetadata]:
    """
    Yield Sentinel-2 tiles that intersect with the given geometries, along with
    their corresponding projection EPSG code and orbit numbers.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    # TODO: some Sentinel-2 tiles slightly overlap UTM zone boundaries, because
    # they are projected and the projected edge of the tile overlaps the UTM
    # zone edge meridian. When the input geometries span multiple UTM zones,
    # *but* also fit entirely within Sentinel-2 projected to a single UTM zone,
    # we should use that single UTM zone.

    # First, find the UTM zones that intersect with the given geometries:
    utm_zones = find_utm_zones(geometries)
    for utm_zone in utm_zones:
        area_of_use = utm_zone.area_of_use
        assert area_of_use

        # Find Sentinel-2 tiles and orbits that intersect with the UTM zone:
        tiles = geopandas.read_file(
            TILES_CONUS_PATH, bbox=area_of_use.bounds, columns=["tile_id"]
        ).set_index("tile_id")
        orbits = geopandas.read_file(
            ORBITS_CONUS_PATH,
            bbox=area_of_use.bounds,
            columns=["relative_orbit_number"],
        ).set_index("relative_orbit_number")

        # Sentinel-2 tiles overlap along UTM zone boundaries. We can narrow
        # down the search by looking at the tile IDs, which start with the UTM
        # zone name:
        zone_name = utm_zone.code.removeprefix("326")
        assert re.match(r"\d{2}", zone_name)
        tiles_in_utm_zone = tiles[tiles.index.str.startswith(zone_name)]

        # Project the tiles and orbits to the UTM zone's CRS:
        tiles_projected = tiles_in_utm_zone.to_crs(epsg=utm_zone.code)
        orbits_projected = orbits.to_crs(epsg=utm_zone.code)

        # The orbit tracks are a line directly beneath the satellite's path.
        # Buffer these tracks to the width of the satellite's field of view:
        orbits_projected.geometry = orbits_projected.geometry.buffer(
            (SWATH_WIDTH_KM + 10) * 1000 / 2,  # add extra 10km to be safe
            cap_style="square",
        )

        # Clip input geometries to UTM zone and project to UTM zone's CRS:
        geometries_clipped_to_utm_zone = geometries.geometry.to_crs("EPSG:4326").clip(
            area_of_use.bounds
        )
        geometries_projected = geometries_clipped_to_utm_zone.to_crs(
            epsg=utm_zone.code
        ).union_all()

        # Select the minimum number of tiles that cover the input geometries:
        tiles_selected = _select_tiles(tiles_projected, geometries_projected)

        orbits_intersection = orbits_projected.geometry.intersection(
            geometries_projected
        )
        orbits_intersecting = orbits_projected.assign(geometry=orbits_intersection)[
            ~orbits_intersection.is_empty
        ]
        for tile_id, relative_orbit_number in (
            tiles_selected.sjoin(orbits_intersecting)
            .drop(columns=["geometry"])
            .itertuples()
        ):
            yield TileMetadata(tile_id, relative_orbit_number)


def find_utm_zones(
    geometries: Union[str, geopandas.GeoDataFrame, geopandas.GeoSeries]
) -> Sequence[CRSInfo]:
    """
    Sentinel-2 rasters are projected using the Universal Transverse Mercator
    (UTM) system.

    Return the UTM zones that intersect with the given geometries.
    """
    if isinstance(geometries, str):
        geometries = geopandas.read_file(geometries)

    geometries = geometries.to_crs("EPSG:4326")
    bounds = bounds_snapped_to_grid(geometries, base=6)

    # Buffer bounds to avoid pulling in adjacent UTM zones:
    bounds["minx"] += 1
    bounds["maxx"] -= 1

    utm_zones = chain.from_iterable(
        query_utm_crs_info(
            datum_name="WGS 84", area_of_interest=AreaOfInterest(minx, miny, maxx, maxy)
        )
        for minx, miny, maxx, maxy in bounds.itertuples(index=False)
    )

    return sorted(
        set(utm_zones),  # deduplicate
        key=lambda zone: int(zone.code),
    )


def _select_tiles(
    tiles: geopandas.GeoDataFrame, geometries_union: shapely.Polygon
) -> geopandas.GeoDataFrame:
    # Find tiles that intersect with the input geometries:
    tiles_intersecting = tiles[tiles.intersects(geometries_union)]

    # Adjacent tiles have a 20km wide overlap where they meet. If the
    # input geometries fall within this boundary area, we don't need to
    # download both tiles. Pick the tiles that cover the most input
    # geometries to minimize the number of rasters we need to download.
    for tile_count in range(1, len(tiles_intersecting) + 1):
        for tile_ids in combinations(tiles_intersecting.index, tile_count):
            tiles = tiles_intersecting.loc[list(tile_ids)]
            if tiles.union_all().covers(geometries_union):
                return tiles
