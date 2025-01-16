import os

import boto3
from botocore import UNSIGNED
from botocore.client import Config

from demeter.raster import Raster
from demeter.raster.usgs.constants import CACHED_RASTER_FILES_DIRECTORY, S3_BUCKET_NAME
from demeter.raster.utils.mask import mask
from demeter.raster.utils.merge import check_for_overlapping_pixels, merge

# Bucket is public, don't send credentials:
s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))


def download_from_s3(key: str) -> str:
    local_path = os.path.join(CACHED_RASTER_FILES_DIRECTORY, key)
    if os.path.exists(local_path):
        # TODO: check if file in cache is stale
        print(f"Cache hit: {local_path}")
    else:
        print(f"Downloading s3://{S3_BUCKET_NAME}/{key}")
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        s3_client.download_file(S3_BUCKET_NAME, key, local_path)

    return local_path


def merge_and_crop_rasters(sources, crop_to=None) -> Raster:
    if crop_to is None:
        return _merge_rasters(sources)

    merged = _merge_rasters(sources, bounds=tuple(crop_to.total_bounds))
    return mask(merged, crop_to, all_touched=True)


def _merge_rasters(sources, **kwargs) -> Raster:
    # TODO: Check that dataset grids line up perfectly. If not, check that
    # merging doesn't cause resampling artifacts.
    # FIXME: merging datasets that are very far apart uses a huge amount of
    # memory, even if the data is very sparse. I think this is because numpy
    # arrays allocate memory for every pixel. Find a way to mitigate.
    print("Merging rasters")

    # USGS elevation tiles overlap their neighboring tiles by 6 pixels. Most of
    # the overlapping data is the same, but neighboring tiles sometimes have
    # different data in the overlapping portion of the boundary.
    # TODO: figure out how to handle this.
    return merge(sources, method=check_for_overlapping_pixels, **kwargs)
