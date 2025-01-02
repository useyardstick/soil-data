import os
import re
from collections.abc import Collection, Iterable, Sequence
from typing import NamedTuple

from defusedxml import ElementTree

from soil_data.raster.sentinel2.constants import Band, Resolution
from soil_data.raster.sentinel2.utils.download import download_keys
from soil_data.raster.sentinel2.utils.search import find_safe_files
from soil_data.raster.sentinel2.utils.tiles import TileMetadata

# Example: "S2B_MSIL2A_20240901T172859_N0511_R055_T14TMM_20240901T215725.SAFE"
SAFE_FILENAME_PATTERN = re.compile(
    r"""\b
    (?P<mission>S2A|S2B)_                   # S2B
    (?P<product_level>MSIL2A)_              # MSIL2A
    (?P<datatake_timestamp>\d{8}T\d{6})_    # 20240901T172859
    (?P<processing_baseline>N\d{4})_        # N0511
    (?P<relative_orbit_number>R\d{3})_      # R055
    (?P<tile_id>T\d{2}[A-Z]{3})_            # T14TMM
    (?P<product_discriminator>\d{8}T\d{6})  # 20240901T215725
    \.SAFE
    \b""",
    re.VERBOSE,
)

# Example: "T14TMM_20240901T172859_B02_10m.jp2"
RASTER_FILENAME_PATTERN = re.compile(
    r"""\b
    (?P<tile_id>T\d{2}[A-Z]{3})_          # T14TMM
    (?P<datatake_timestamp>\d{8}T\d{6})_  # 20240901T172859
    (?P<band>[A-Z\d]{3})_                 # B02
    (?P<resolution>\d+?m)                 # 10m
    \.jp2
    $""",
    re.VERBOSE,
)

# Example: "MSK_DETFOO_B02.jp2"
DETECTOR_FOOTPRINT_MASK_FILENAME_PATTERN = re.compile(
    r"\bMSK_DETFOO_(?P<band>[A-Z\d]{3})\.jp2$"
)


class SafeMetadata(NamedTuple):
    tile_id: str
    datatake_timestamp: str

    @classmethod
    def from_filename(cls, name: str) -> "SafeMetadata":
        match = re.search(SAFE_FILENAME_PATTERN, name)
        if match is None:
            raise ValueError(f"Could not parse SAFE metadata from filename: {name}")

        return cls(
            tile_id=match.group("tile_id").removeprefix("T"),
            datatake_timestamp=match.group("datatake_timestamp"),
        )

    @property
    def utm_zone(self) -> str:
        zone_name = self.tile_id[:2]
        assert re.match(r"\d{2}", zone_name)
        return zone_name

    @property
    def crs(self) -> str:
        return f"EPSG:326{self.utm_zone}"


class RasterMetadata(NamedTuple):
    tile_id: str
    datatake_timestamp: str
    band: str
    resolution: str

    @classmethod
    def from_filename(cls, name: str) -> "RasterMetadata":
        match = re.search(RASTER_FILENAME_PATTERN, name)
        if match is None:
            raise ValueError(f"Could not parse raster metadata from filename: {name}")

        return cls(
            tile_id=match.group("tile_id").removeprefix("T"),
            datatake_timestamp=match.group("datatake_timestamp"),
            band=match.group("band"),
            resolution=match.group("resolution"),
        )


class DetectorFootprintMaskMetadata(NamedTuple):
    band: str

    @classmethod
    def from_filename(cls, name: str) -> "DetectorFootprintMaskMetadata":
        match = re.search(DETECTOR_FOOTPRINT_MASK_FILENAME_PATTERN, name)
        if match is None:
            raise ValueError(
                f"Could not parse detector footprint mask metadata from filename: {name}"
            )

        return cls(band=match.group("band"))


def list_raster_keys(
    safe_keys: Sequence[str],
    bands: Collection[tuple[Band, Resolution]],
) -> Iterable[str]:
    """
    Search the given SAFE file manifests for rasters for specific spectral
    bands. Yield the S3 key for each raster file found.
    """
    bands_to_download = {(band.value, resolution.value) for band, resolution in bands}
    band_masks_to_download = {band.value for band, _ in bands}

    # Download the manifest file for each tile:
    manifest_paths = download_keys(f"{key}/manifest.safe" for key in safe_keys)

    # Parse each manifest's xml:
    for safe_key, manifest_path in zip(safe_keys, manifest_paths):
        for raster_key in extract_raster_keys_from_manifest(safe_key, manifest_path):
            # For each file path, first check if it's matches the raster
            # filename pattern:
            try:
                raster_metadata = RasterMetadata.from_filename(raster_key)
            except ValueError:
                # If not, try the detector footprint mask filename pattern:
                try:
                    mask_metadata = DetectorFootprintMaskMetadata.from_filename(
                        raster_key
                    )
                except ValueError:
                    # If it's not a raster or mask file, skip it
                    continue
                # If it's a detector footprint mask, check if it's for one of
                # the requested bands:
                if mask_metadata.band in band_masks_to_download:
                    yield raster_key
            else:
                # If it's a data raster, check if it's for one of the requested
                # bands, at the requested resolution:
                if (
                    raster_metadata.band,
                    raster_metadata.resolution,
                ) in bands_to_download:
                    yield raster_key


def extract_raster_keys_from_manifest(
    safe_key: str, manifest_path: str
) -> Iterable[str]:
    manifest_tree = ElementTree.parse(manifest_path)

    # Scan every file path listed in the manifest:
    for element in manifest_tree.iter("fileLocation"):
        raster_path = element.get("href")
        if raster_path is None:
            continue

        yield os.path.normpath(os.path.join(safe_key, raster_path))
