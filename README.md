Demeter
=======

Tools for fetching spatial datasets for agriculture.

Usage
-----

See API documentation:

- [POLARIS](docs.md#demeter.raster.polaris)
- [Soil and Landscape Grid of Australia (SLGA)](docs.md#demeterrasterslga)
- USGS
  - [Topography](docs.md#demeterrasterusgstopography)
  - [Hydrography](docs.md#demeterrasterusgshydrography)
- [NDVI from Sentinel-2 imagery](docs.md##demeterrastersentinel2ndvi)

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
