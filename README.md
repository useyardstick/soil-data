Demeter
=======

Tools for fetching spatial datasets for agriculture.

Installation
------------

```
pip install git+https://github.com/useyardstick/demeter.git
```

Usage
-----

### Fetching raster datasets

- [POLARIS](docs.md#demeter.raster.polaris)
- [Soil and Landscape Grid of Australia (SLGA)](docs.md#demeter.raster.slga)
- USGS
  - [Topography](docs.md#demeter.raster.usgs.topography)
  - [Hydrography](docs.md#demeter.raster.usgs.hydrography)
- [NDVI from Sentinel-2 imagery](docs.md##demeter.raster.sentinel2.ndvi)

Tools in this library return [`Raster`](docs.md#demeter.raster.Raster)
instances. To save to disk, use the `save` method:

```python
raster.save("path/to/file.tif")
```

#### Raster helpers

This library also provides helper functions for common raster operations:

- [Masking](docs.md#demeter.raster.utils.mask)
- [Merging](docs.md#demeter.raster.utils.merge)

Development
-----------

Create a virtual environment (recommended):

```
python -m venv .venv
source .venv/bin/activate
```

Install development dependencies:

```
pip install -r requirements.dev.txt
```

Install pre-commit hooks:

```
pre-commit install
```

Run tests:

```
pytest
```
