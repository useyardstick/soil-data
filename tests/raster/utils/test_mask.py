import geopandas
import numpy
from rasterio import Affine

from demeter.raster import Raster
from demeter.raster.utils.mask import mask_raster


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
        Raster(matrix, transform=Affine.identity(), crs="EPSG:4326"),
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
