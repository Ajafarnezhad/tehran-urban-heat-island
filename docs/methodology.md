# Scientific Methodology

## 1. Study Area
The Tehran metropolitan area (approximate extent: `51.05°E–51.60°E`, `35.55°N–35.85°N`) was chosen as one of the most densely built-up metropolitan areas in the Middle East, characterized by rapid urban expansion and shrinking green space.

## 2. Data

| Source | Role | Spatial resolution | Time period |
|---|---|---|---|
| Landsat 8 Collection 2 Level-2 | LST/NDVI/NDBI retrieval, baseline period | 30 m (thermal band resampled) | Summer 2015 |
| Landsat 9 Collection 2 Level-2 | LST/NDVI/NDBI retrieval, recent period | 30 m | Summer 2024 |
| ESA WorldCover v200 | Land-use / land-cover classification | 10 m | 2021 |

A cloud-cover filter below 20% and a cloud/shadow mask from the `QA_PIXEL` band (CFMask bits) are applied to ensure thermal pixel quality.

## 3. Land Surface Temperature (LST) Retrieval

1. Convert raw thermal band (`ST_B10`) values to Kelvin using the official USGS Collection 2 calibration coefficients.
2. Estimate land surface emissivity (LSE) using the **NDVI Threshold Method**:
   - Water (`NDVI<0`): ε = 0.991
   - Bare soil (`0≤NDVI<0.2`): ε = 0.966
   - Full vegetation cover (`NDVI>0.5`): ε = 0.973
   - Mixed pixel: interpolated from the vegetation fraction (Pv), plus a surface-roughness correction term (0.005)
3. Correct for emissivity using the simplified Planck equation (Artis & Carnahan, 1982) to convert brightness temperature into true LST.
4. Convert the final result from Kelvin to Celsius.

## 4. Complementary Spectral Indices
- **NDVI** (greenness): `(NIR−Red)/(NIR+Red)`
- **NDBI** (built-up density): `(SWIR1−NIR)/(SWIR1+NIR)`
- **NDWI** (surface water, used to exclude water pixels from the regression analysis): `(Green−NIR)/(Green+NIR)`

## 5. Statistical and Spatial Analysis

- **Grid sampling**: a regular 500×500 m grid is generated over the study area, and the LST/NDVI/NDBI/NDWI bands are zonally averaged within each cell.
- **Pearson correlation**: measures the strength and direction of the relationship between LST and the spectral indices, quantitatively demonstrating the "vegetation cooling effect."
- **Getis-Ord Gi\* statistic**: identifies statistically significant clusters of high LST (hot spots) and low LST (cold spots) using a Queen contiguity weights matrix and 999 random permutations; classified at three confidence levels (90%, 95%, 99%).
- **UHI intensity**: the difference between the mean LST of "built-up" areas (WorldCover class 50) and the mean LST of natural/agricultural land cover (classes 10, 20, 30, 40) is used as a quantitative measure of heat island intensity.
- **Temporal trend analysis**: LST distributions and means are compared between summer 2015 and summer 2024 to estimate the urban warming rate over the nine-year interval.

## 6. Real-Data Validation Run

In addition to the general Earth Engine pipeline (any city, any date range), this project includes a fully reproducible run against two specific, real Landsat scenes, fetched directly from the public [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) STAC catalog — no Google account or cloud project required:

| Scene | Date | Cloud cover |
|---|---|---|
| `LC08_L2SP_164035_20150805_02_T1` | 2015-08-05 | 0.93% |
| `LC09_L2SP_164035_20240805_02_T1` | 2024-08-05 | 0.53% |

Both scenes were deliberately picked on the same calendar day, nine years apart, to control for seasonal sun-angle effects as much as a two-image comparison allows. Results are in [README.md](../README.md#real-results).

## 7. Limitations
- The mono-window method is somewhat less accurate than multi-channel (split-window) methods, but is simpler to apply with a single thermal band, as on Landsat.
- Differences in satellite overpass time and atmospheric conditions between the two periods may explain part of the observed temperature difference; averaging several scenes within a similar seasonal window helps but does not eliminate this. In the real-data validation run above, a single-date-pair comparison showed a negligible 0.05°C difference in mean LST — confirming that daily weather noise dominates over any long-term warming signal at this sample size, and that a robust trend estimate would require multi-year compositing, not a two-scene snapshot.
- 30 m spatial resolution for LST is not sufficient for precise intra-neighborhood analysis; higher-resolution thermal data (such as ECOSTRESS, ~70 m) could reduce this limitation in future iterations.

## 8. Key References
- Artis, D.A. & Carnahan, W.H. (1982). Survey of emissivity variability in thermography of urban areas. *Remote Sensing of Environment*.
- Sobrino, J.A. et al. (2004). Land surface temperature retrieval from LANDSAT TM 5. *Remote Sensing of Environment*.
- Getis, A. & Ord, J.K. (1992). The analysis of spatial association by use of distance statistics. *Geographical Analysis*.
- Zhou, D. et al. (2018). Satellite remote sensing of surface urban heat islands: progress, challenges, and perspectives. *Remote Sensing*.
