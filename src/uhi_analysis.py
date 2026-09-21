"""Statistical and spatial analysis of the heat island: grid sampling, correlation, and Getis-Ord Gi* hotspots."""
import ee
import numpy as np
import pandas as pd
import geopandas as gpd
from libpysal.weights import Queen
from esda.getisord import G_Local


def sample_grid_to_dataframe(
    composite: "ee.Image",
    aoi: "ee.Geometry",
    cell_size_m: int = 500,
    scale: int = 30,
) -> pd.DataFrame:
    """Build a regular grid over the AOI and extract the mean LST/NDVI/NDBI/NDWI per cell."""
    grid = composite.select("LST_C").reduceRegions(
        collection=ee.FeatureCollection([ee.Feature(aoi)]).geometry().coveringGrid(
            ee.Projection("EPSG:32639"), cell_size_m
        ),
        reducer=ee.Reducer.mean(),
        scale=scale,
    )

    full = composite.reduceRegions(
        collection=grid,
        reducer=ee.Reducer.mean(),
        scale=scale,
    )

    features = full.getInfo()["features"]
    rows = []
    for feat in features:
        props = feat["properties"]
        geom = feat["geometry"]
        rows.append({**props, "geometry": geom})

    df = pd.DataFrame(rows)
    return df


def to_geodataframe(df: pd.DataFrame, crs: str = "EPSG:4326") -> gpd.GeoDataFrame:
    from shapely.geometry import shape

    geometries = [shape(g) for g in df["geometry"]]
    gdf = gpd.GeoDataFrame(df.drop(columns=["geometry"]), geometry=geometries, crs=crs)
    return gdf.dropna(subset=["LST_C"])


def correlation_report(gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Pearson correlation matrix between LST and spectral indices, quantifying the greening-vs-heat relationship."""
    cols = [c for c in ["LST_C", "NDVI", "NDBI", "NDWI"] if c in gdf.columns]
    return gdf[cols].corr(method="pearson")


def getis_ord_hotspots(gdf: gpd.GeoDataFrame, value_col: str = "LST_C") -> gpd.GeoDataFrame:
    """Identify statistically significant hot spot / cold spot clusters using the Getis-Ord Gi* statistic.

    The output includes a z-score column and a confidence classification (90/95/99%) for mapping.
    """
    gdf = gdf.reset_index(drop=True)
    w = Queen.from_dataframe(gdf, use_index=False)
    w.transform = "r"

    y = gdf[value_col].values
    g_local = G_Local(y, w, star=True, permutations=999)

    gdf["gi_zscore"] = g_local.Zs
    gdf["gi_pvalue"] = g_local.p_sim

    def classify(z, p):
        if p > 0.10:
            return "Not significant"
        if z > 0:
            if p <= 0.01:
                return "Hot spot - 99% confidence"
            if p <= 0.05:
                return "Hot spot - 95% confidence"
            return "Hot spot - 90% confidence"
        else:
            if p <= 0.01:
                return "Cold spot - 99% confidence"
            if p <= 0.05:
                return "Cold spot - 95% confidence"
            return "Cold spot - 90% confidence"

    gdf["hotspot_class"] = [classify(z, p) for z, p in zip(gdf["gi_zscore"], gdf["gi_pvalue"])]
    return gdf


def uhi_intensity(gdf: gpd.GeoDataFrame, worldcover_col: str = "worldcover") -> dict:
    """Compute UHI intensity: the difference in mean LST between built-up areas and a green/natural reference."""
    urban_mean = gdf.loc[gdf[worldcover_col] == 50, "LST_C"].mean()
    rural_mean = gdf.loc[gdf[worldcover_col].isin([10, 20, 30, 40]), "LST_C"].mean()
    return {
        "urban_mean_lst_c": round(float(urban_mean), 2),
        "rural_mean_lst_c": round(float(rural_mean), 2),
        "uhi_intensity_c": round(float(urban_mean - rural_mean), 2),
    }
