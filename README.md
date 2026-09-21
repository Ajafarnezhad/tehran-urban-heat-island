# 🛰️ Tehran Urban Heat Island Analysis (Remote Sensing)

A complete, reproducible pipeline for retrieving Land Surface Temperature (LST) from Landsat satellite imagery, evaluating its relationship with vegetation and built-up surfaces, identifying statistically significant heat clusters, and comparing Tehran's urban warming pattern between 2015 and 2024.

Built on **Google Earth Engine (Python API)**, **geemap**, and standard spatial data science libraries (`geopandas`, `esda`, `libpysal`) — with a second, no-account-required execution path against real satellite imagery.

---

## Why this project

The urban heat island effect is one of the most tangible consequences of unchecked urban growth and vanishing green space. This project shows how, using nothing but free satellite imagery and no field data, that effect can be measured quantitatively, analyzed spatially, and visualized clearly — a complete demonstration of the remote sensing pipeline from raw imagery to spatial statistics and urban decision support.

## Key features

- ✅ LST retrieval via the mono-window method with NDVI-based emissivity correction
- ✅ Fully automated CLI pipeline, runnable with a single command
- ✅ Two independent execution paths: Google Earth Engine (any city/date range) or a no-account fallback using real, pre-selected Landsat scenes
- ✅ Temporal comparison (2015 vs. 2024) to estimate the urban warming rate
- ✅ Correlation analysis between LST and NDVI/NDBI/NDWI
- ✅ Statistically significant hot/cold spot detection with **Getis-Ord Gi\*** (genuine spatial statistics, not just visual inspection)
- ✅ Quantitative UHI intensity metric (temperature gap between built-up and natural areas)
- ✅ Interactive web map (Leaflet/geemap) and report-ready charts
- ✅ Modular code, portable to any other city by editing a single configuration file

## Project layout

```
heat/
├── config/config.yaml           # AOI, time periods, analysis parameters
├── src/
│   ├── gee_setup.py              # Earth Engine connection and authentication
│   ├── indices.py                 # Cloud masking, NDVI, NDBI, NDWI
│   ├── lst.py                     # LST retrieval (emissivity + mono-window)
│   ├── lulc.py                    # ESA WorldCover land cover
│   ├── uhi_analysis.py            # Grid sampling, correlation, Getis-Ord Gi*
│   ├── visualization.py           # Interactive maps and charts
│   └── export.py                  # GeoTIFF/CSV export to Google Drive
├── scripts/
│   ├── run_pipeline.py           # Path 1: full run on Earth Engine
│   ├── fetch_real_scenes.py      # Path 2: download real scenes (no login)
│   └── analyze_real_scenes.py    # Path 2: analyze that real data
├── notebooks/01_UHI_Tehran_Analysis.ipynb   # Step-by-step analytical narrative
├── docs/methodology.md           # Full scientific methodology + references
├── docs/portfolio_page.html      # Standalone case-study page with real results
└── outputs/                      # Maps, charts, and result tables
```

## Setup

```bash
pip install -r requirements.txt
```

## Running it — two paths

**Path 1: Google Earth Engine (any city, any time range)**

Create a Google Cloud project with the Earth Engine API enabled (free for non-commercial use) and put its ID in `config/config.yaml` (`project.ee_project_id`). The first run opens a browser window for Google sign-in.

```bash
python scripts/run_pipeline.py
# or, for the interactive narrative:
jupyter lab notebooks/01_UHI_Tehran_Analysis.ipynb
```

**Path 2: Reproduce without any account (real data bundled with this repo)**

This path downloads two real Landsat scenes directly from the public [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) catalog — no Google authentication needed:

```bash
python scripts/fetch_real_scenes.py     # download real bands over Tehran
python scripts/analyze_real_scenes.py   # run the same LST/NDVI/Getis-Ord algorithm on real data
```

## Real results (produced by Path 2, on real Landsat 8/9 scenes over Tehran)

| Scene | Date | Cloud cover |
|---|---|---|
| `LC08_L2SP_164035_20150805_02_T1` | 2015-08-05 | 0.93% |
| `LC09_L2SP_164035_20240805_02_T1` | 2024-08-05 | 0.53% |

Across 2,666 valid 500 m grid cells:

| Metric | Value |
|---|---|
| Mean LST, low-vegetation areas (NDVI<0.15) | 50.34 °C |
| Mean LST, green areas (NDVI>0.35) | 45.42 °C |
| Urban Heat Island intensity | **4.92 °C** |
| Pearson correlation, NDVI–LST | r = −0.295 (p < 10⁻⁵⁰) |
| Pearson correlation, NDBI–LST | r = 0.552 (p < 10⁻²⁰⁰) |
| Mean LST difference between the two dates (single-day snapshot) | +0.05 °C — not significant; a robust long-term trend needs multi-year compositing |

A Getis-Ord Gi* analysis on the same grid classified 830 cells as "hot spot" and 1,018 cells as "cold spot" at ≥90% confidence (full detail: `outputs/data/summary.json`).

A visual version of these results (charts and maps) is available in the [project case-study page](docs/portfolio_page.html).

## Sample outputs (in `outputs/` after running)

| File | Description |
|---|---|
| `figures/lst_histogram.png` | Temperature distribution, period comparison |
| `figures/ndvi_lst_scatter.png` | Vegetation vs. temperature relationship |
| `figures/hotspot_map.png` | Getis-Ord Gi* hot/cold spot map |
| `figures/correlation_heatmap.png` | Index correlation matrix |
| `data/summary.json` | Numeric summary of results (mean temperatures, UHI intensity, correlations) |
| `data/hotspots.geojson` | Hotspot spatial layer for GIS |

## Adapting to other cities

Just change `bbox` in `config/config.yaml` to the target city's extent — the entire pipeline runs without any code changes.

## Scientific basis

Full methodology, formulas, and references are in [`docs/methodology.md`](docs/methodology.md).

## Tech stack

`Google Earth Engine` · `Microsoft Planetary Computer / STAC` · `Python` · `geemap` · `rasterio` · `geopandas` · `esda / libpysal (PySAL)` · `scipy` · `matplotlib / seaborn` · `Landsat 8/9 Collection 2` · `ESA WorldCover`

## Author

**Amirhossein Jafarnezhad**

## License

MIT — free for educational and research use.
