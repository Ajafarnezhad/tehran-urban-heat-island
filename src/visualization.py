"""ساخت نقشه‌های تعاملی (geemap/folium) و نمودارهای ثابت (matplotlib) برای گزارش."""
import geemap
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
import pandas as pd

LST_VIS = {"min": 20, "max": 55, "palette": ["040274", "3ae237", "ffd611", "ff8b13", "ff0000", "990000"]}
NDVI_VIS = {"min": -0.2, "max": 0.8, "palette": ["a50026", "ffffbf", "1a9850"]}


def make_interactive_map(aoi, lst_baseline, lst_recent, worldcover, center_zoom=11):
    m = geemap.Map()
    m.centerObject(aoi, center_zoom)
    m.addLayer(lst_baseline.select("LST_C"), LST_VIS, "LST 2015")
    m.addLayer(lst_recent.select("LST_C"), LST_VIS, "LST 2024")
    m.addLayer(lst_recent.select("NDVI"), NDVI_VIS, "NDVI 2024")
    m.addLayer(worldcover.randomVisualizer(), {}, "پوشش اراضی (ESA WorldCover)")
    m.add_colorbar(LST_VIS, label="دمای سطح زمین (°C)")
    m.addLayerControl()
    return m


def plot_lst_histogram(gdf_baseline: pd.DataFrame, gdf_recent: pd.DataFrame, out_path: str):
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(gdf_baseline["LST_C"], label="تابستان 2015", fill=True, ax=ax)
    sns.kdeplot(gdf_recent["LST_C"], label="تابستان 2024", fill=True, ax=ax)
    ax.set_xlabel("دمای سطح زمین (°C)")
    ax.set_ylabel("چگالی")
    ax.set_title("توزیع دمای سطح زمین در تهران — مقایسه دو دوره")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_ndvi_lst_scatter(gdf: pd.DataFrame, out_path: str):
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.regplot(data=gdf, x="NDVI", y="LST_C", scatter_kws={"alpha": 0.4, "s": 12}, line_kws={"color": "red"}, ax=ax)
    ax.set_title("رابطه بین پوشش گیاهی (NDVI) و دمای سطح زمین")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_hotspot_map(gdf: gpd.GeoDataFrame, out_path: str):
    fig, ax = plt.subplots(figsize=(8, 8))
    gdf.plot(column="hotspot_class", categorical=True, legend=True, ax=ax, cmap="RdYlBu_r", edgecolor="white", linewidth=0.1)
    ax.set_title("خوشه‌های داغ و سرد دمایی تهران (آماره Getis-Ord Gi*)")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_correlation_heatmap(corr_df: pd.DataFrame, out_path: str):
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(corr_df, annot=True, cmap="coolwarm", vmin=-1, vmax=1, ax=ax)
    ax.set_title("ماتریس همبستگی شاخص‌ها")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
