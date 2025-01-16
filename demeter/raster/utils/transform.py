from rasterio import Affine


def extract_resolution_from_transform(
    transform: Affine,
) -> tuple[float, float]:
    """
    Return (x, y) resolution of the given Affine transform.
    """
    return transform.a, -transform.e
