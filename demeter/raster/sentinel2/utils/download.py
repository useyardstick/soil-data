import os
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import boto3

from demeter.raster.sentinel2.constants import S3_BUCKET_NAME

S3_PREFIX = "Sentinel-2/MSI/L2A/"


def get_cache_directory():
    """
    Downloaded raster files are cached here.
    """
    return os.environ.get("SENTINEL2_CACHED_RASTER_FILES_DIRECTORY", ".sentinel2_cache")


def download_keys(keys: Iterable[str]) -> Iterable[str]:
    """
    Download all the Sentinel-2 rasters with the given keys. Yield the local
    path for each downloaded file.
    """
    s3_client = _s3_client()
    cache_directory = get_cache_directory()
    downloader = partial(_download_from_s3, s3_client, cache_directory)

    # We don't really need concurrency here, since boto3 already downloads
    # large files in parallel chunks. This ThreadPoolExecutor is just a
    # convenient way to move downloading to the background.
    pool = ThreadPoolExecutor(max_workers=1)

    try:
        results = pool.map(downloader, keys)
    finally:
        # Don't wait until all downloads are complete before returning the
        # results. This allows callers to start processing the results while
        # the downloads are still in progress:
        pool.shutdown(wait=False)

    return results


def _s3_client():
    """
    Client for Copernicus' S3 API. Note that this is not AWS S3.

    Manage credentials here (Chris has login in 1Password):
    https://eodata-s3keysmanager.dataspace.copernicus.eu/panel/s3-credentials
    """
    return boto3.client(
        "s3",
        endpoint_url=os.environ.get("COPERNICUS_AWS_ENDPOINT_URL"),
        aws_access_key_id=os.environ["COPERNICUS_AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["COPERNICUS_AWS_SECRET_ACCESS_KEY"],
    )


def _download_from_s3(s3_client, cache_directory, key: str) -> str:
    key = os.path.normpath(key)
    assert key.startswith(S3_PREFIX)

    local_path = os.path.join(cache_directory, key)
    if os.path.exists(local_path):
        # TODO: verify checksum from manifest to see if file in cache is stale
        print(f"Cache hit: {local_path}")
    else:
        print(f"Downloading s3://{S3_BUCKET_NAME}/{key}")
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        s3_client.download_file(S3_BUCKET_NAME, key, local_path)

    return local_path
