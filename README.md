# Azimuthal anisotropy from inter-station velocity measurements

[![DOI](https://zenodo.org/badge/264389474.svg)](https://zenodo.org/badge/latestdoi/264389474)

Determine the 2-Theta and 4-Theta terms of azimuthal anisotropy from inter-station group-velocity residuals in overlapping cells. This methodology is established in Schippkus et al. (2020), but the original version used in the paper requires the specific output of a specific isotropic inversion code, relying on the parametrization of the isotropic inversion. This new version is based purely on a .csv file of station-pair locations and inter-station group-velocity residuals.

Format of the header-less `.csv` input: `lat1, lon1, lat2, lon2, vel`

[1] Schippkus, S., Zigone, D., Bokelmann, G., AlpArray Working Group. (2020). Azimuthal anisotropy in the wider Vienna basin region: a proxy for the present-day stress field and deformation. Geophysical Journal International, 220(3), 2056–2067. http://doi.org/10.1093/gji/ggz565

## Requirements

- Python 3.5+
- Numpy
- Matplotlib
- Obspy (https://github.com/obspy/obspy/wiki)
- Shapely (https://pypi.org/project/Shapely/)

## TODO

- [x] Implement weighted average using path-lengths within cells.
- [ ] Actually test with real data.
- [ ] Synthetic example that can handle varying anisotropy.
