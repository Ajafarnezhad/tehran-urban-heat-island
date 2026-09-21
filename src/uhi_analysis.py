"""تحلیل آماری و فضایی جزیره حرارتی: نمونه‌برداری شبکه‌ای، همبستگی و هات‌اسپات Getis-Ord Gi*."""
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
    """تولید شبکه منظم روی AOI و استخراج میانگین LST/NDVI/NDBI/NDWI در هر سلول."""
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
    """ماتریس همبستگی پیرسون بین LST و شاخص‌های طیفی، برای گزارش کمّی رابطه سبزینگی-حرارت."""
    cols = [c for c in ["LST_C", "NDVI", "NDBI", "NDWI"] if c in gdf.columns]
    return gdf[cols].corr(method="pearson")


def getis_ord_hotspots(gdf: gpd.GeoDataFrame, value_col: str = "LST_C") -> gpd.GeoDataFrame:
    """شناسایی خوشه‌های داغ (Hot Spot) و سرد (Cold Spot) آماری با آماره Getis-Ord Gi*.

    خروجی شامل ستون z_score و طبقه‌بندی اطمینان (90/95/99٪) برای نقشه‌سازی است.
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
            return "بدون معناداری آماری"
        if z > 0:
            if p <= 0.01:
                return "کانون داغ - اطمینان 99٪"
            if p <= 0.05:
                return "کانون داغ - اطمینان 95٪"
            return "کانون داغ - اطمینان 90٪"
        else:
            if p <= 0.01:
                return "کانون سرد - اطمینان 99٪"
            if p <= 0.05:
                return "کانون سرد - اطمینان 95٪"
            return "کانون سرد - اطمینان 90٪"

    gdf["hotspot_class"] = [classify(z, p) for z, p in zip(gdf["gi_zscore"], gdf["gi_pvalue"])]
    return gdf


def uhi_intensity(gdf: gpd.GeoDataFrame, worldcover_col: str = "worldcover") -> dict:
    """محاسبه شدت جزیره حرارتی: اختلاف میانگین LST مناطق ساخته‌شده و مناطق سبز/طبیعی مرجع."""
    urban_mean = gdf.loc[gdf[worldcover_col] == 50, "LST_C"].mean()
    rural_mean = gdf.loc[gdf[worldcover_col].isin([10, 20, 30, 40]), "LST_C"].mean()
    return {
        "urban_mean_lst_c": round(float(urban_mean), 2),
        "rural_mean_lst_c": round(float(rural_mean), 2),
        "uhi_intensity_c": round(float(urban_mean - rural_mean), 2),
    }
