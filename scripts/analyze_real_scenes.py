"""Real analysis of the Tehran urban heat island on real downloaded Landsat bands.

Input: outputs/data/raw_bands/{baseline,recent}.npz (produced by fetch_real_scenes.py)
Implements the LST algorithm from src/lst.py (mono-window method + NDVI emissivity) using
numpy instead of Earth Engine, since this script runs on pre-downloaded data rather than
the cloud service.

Outputs are saved to outputs/figures and outputs/data and are the real source of the numbers
reported in the README and the project's portfolio page.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from shapely.geometry import box
from scipy import stats
from libpysal.weights import Queen
from esda.getisord import G_Local

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "outputs" / "data" / "raw_bands"
FIG_DIR = ROOT / "outputs" / "figures"
DATA_DIR = ROOT / "outputs" / "data"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SCALE_OPTICAL_MULT, SCALE_OPTICAL_ADD = 0.0000275, -0.2
SCALE_THERMAL_MULT, SCALE_THERMAL_ADD = 0.00341802, 149.0

QA_CLOUD_SHADOW_BIT = 1 << 4
QA_CLOUD_BIT = 1 << 3
QA_CIRRUS_BIT = 1 << 2


def load_scene(period_key: str):
    npz = np.load(RAW_DIR / f"{period_key}.npz")
    meta = json.load(open(RAW_DIR / f"{period_key}_meta.json", encoding="utf-8"))
    transform = npz["transform"]
    return {
        "red": npz["red"].astype(np.float64),
        "nir": npz["nir08"].astype(np.float64),
        "swir1": npz["swir16"].astype(np.float64),
        "green": npz["green"].astype(np.float64),
        "thermal": npz["lwir11"].astype(np.float64),
        "qa": npz["qa_pixel"].astype(np.int32),
        "transform": transform,
        "meta": meta,
    }


def cloud_mask(qa: np.ndarray) -> np.ndarray:
    bad = (
        (qa & QA_CLOUD_SHADOW_BIT != 0)
        | (qa & QA_CLOUD_BIT != 0)
        | (qa & QA_CIRRUS_BIT != 0)
    )
    return ~bad


def compute_layers(scene: dict) -> dict:
    valid = cloud_mask(scene["qa"]) & (scene["red"] > 0) & (scene["thermal"] > 0)

    red = scene["red"] * SCALE_OPTICAL_MULT + SCALE_OPTICAL_ADD
    nir = scene["nir"] * SCALE_OPTICAL_MULT + SCALE_OPTICAL_ADD
    swir1 = scene["swir1"] * SCALE_OPTICAL_MULT + SCALE_OPTICAL_ADD
    green = scene["green"] * SCALE_OPTICAL_MULT + SCALE_OPTICAL_ADD
    thermal_k = scene["thermal"] * SCALE_THERMAL_MULT + SCALE_THERMAL_ADD

    with np.errstate(invalid="ignore", divide="ignore"):
        ndvi = (nir - red) / (nir + red)
        ndbi = (swir1 - nir) / (swir1 + nir)
        ndwi = (green - nir) / (green + nir)

        pv = np.clip((ndvi - 0.2) / 0.3, 0, 1) ** 2
        emissivity = np.where(
            ndvi < 0, 0.991,
            np.where(ndvi < 0.2, 0.966,
                     np.where(ndvi > 0.5, 0.973, 0.966 * (1 - pv) + 0.973 * pv + 0.005)),
        )

        wavelength = 10.895
        rho = 14388.0
        lst_k = thermal_k / (1 + (wavelength * thermal_k / rho) * np.log(emissivity))
        lst_c = lst_k - 273.15

    for arr in (ndvi, ndbi, ndwi, lst_c):
        arr[~valid] = np.nan
    lst_c[(lst_c < 0) | (lst_c > 70)] = np.nan

    return {"ndvi": ndvi, "ndbi": ndbi, "ndwi": ndwi, "lst_c": lst_c, "valid": valid}


def grid_to_geodataframe(layers: dict, transform, crs: str, cell_px: int = 17) -> gpd.GeoDataFrame:
    """Aggregate 30 m pixels into ~500 m cells (17 pixels) and build a GeoDataFrame."""
    lst = layers["lst_c"]
    ndvi = layers["ndvi"]
    ndbi = layers["ndbi"]
    h, w = lst.shape

    rows = []
    for i in range(0, h - cell_px, cell_px):
        for j in range(0, w - cell_px, cell_px):
            block_lst = lst[i:i + cell_px, j:j + cell_px]
            block_ndvi = ndvi[i:i + cell_px, j:j + cell_px]
            block_ndbi = ndbi[i:i + cell_px, j:j + cell_px]
            if np.isnan(block_lst).mean() > 0.4:
                continue
            x0, y0 = transform * (j, i)
            x1, y1 = transform * (j + cell_px, i + cell_px)
            rows.append({
                "LST_C": np.nanmean(block_lst),
                "NDVI": np.nanmean(block_ndvi),
                "NDBI": np.nanmean(block_ndbi),
                "geometry": box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)),
            })
    return gpd.GeoDataFrame(rows, crs=crs)


def getis_ord(gdf: gpd.GeoDataFrame, value_col="LST_C") -> gpd.GeoDataFrame:
    gdf = gdf.reset_index(drop=True)
    w = Queen.from_dataframe(gdf, use_index=False)
    w.transform = "r"
    g = G_Local(gdf[value_col].values, w, star=True, permutations=999, seed=42)
    gdf["gi_zscore"] = g.Zs
    gdf["gi_pvalue"] = g.p_sim

    def classify(z, p):
        if p > 0.10:
            return "Not significant"
        conf = "99%" if p <= 0.01 else ("95%" if p <= 0.05 else "90%")
        kind = "Hot spot" if z > 0 else "Cold spot"
        return f"{kind} - {conf} confidence"

    gdf["hotspot_class"] = [classify(z, p) for z, p in zip(gdf["gi_zscore"], gdf["gi_pvalue"])]
    return gdf


def main():
    print("Loading real scenes...")
    baseline_raw = load_scene("baseline")
    recent_raw = load_scene("recent")

    print("Computing NDVI/NDBI/LST from real pixels...")
    baseline = compute_layers(baseline_raw)
    recent = compute_layers(recent_raw)

    transform = recent_raw["transform"]
    from affine import Affine
    affine_t = Affine(*transform)

    print("Building analysis grid (~500 m)...")
    gdf_recent = grid_to_geodataframe(recent, affine_t, recent_raw["meta"]["crs"])
    gdf_baseline = grid_to_geodataframe(baseline, affine_t, baseline_raw["meta"]["crs"])

    gdf_recent = gdf_recent.dropna(subset=["LST_C", "NDVI"])
    gdf_baseline = gdf_baseline.dropna(subset=["LST_C", "NDVI"])

    print(f"Valid cells — baseline: {len(gdf_baseline)}, recent: {len(gdf_recent)}")

    print("Computing Pearson correlations...")
    valid_pairs = gdf_recent[["NDVI", "NDBI", "LST_C"]].dropna()
    r_ndvi, p_ndvi = stats.pearsonr(valid_pairs["NDVI"], valid_pairs["LST_C"])
    r_ndbi, p_ndbi = stats.pearsonr(valid_pairs["NDBI"], valid_pairs["LST_C"])
    corr_matrix = valid_pairs.corr(method="pearson")
    corr_matrix.to_csv(DATA_DIR / "correlation_matrix.csv", encoding="utf-8-sig")

    print("Running Getis-Ord Gi*...")
    gdf_hotspots = getis_ord(gdf_recent, "LST_C")
    gdf_hotspots.to_file(DATA_DIR / "hotspots.geojson", driver="GeoJSON")
    gdf_recent.to_file(DATA_DIR / "grid_recent.geojson", driver="GeoJSON")
    gdf_baseline.to_file(DATA_DIR / "grid_baseline.geojson", driver="GeoJSON")

    ndvi_thresh = 0.35
    urban_mask_recent = gdf_recent["NDVI"] < 0.15
    green_mask_recent = gdf_recent["NDVI"] > ndvi_thresh
    urban_mean = gdf_recent.loc[urban_mask_recent, "LST_C"].mean()
    green_mean = gdf_recent.loc[green_mask_recent, "LST_C"].mean()
    uhi_intensity = urban_mean - green_mean

    warming = gdf_recent["LST_C"].mean() - gdf_baseline["LST_C"].mean()

    summary = {
        "baseline_scene": baseline_raw["meta"]["scene_id"],
        "baseline_date": baseline_raw["meta"]["date"],
        "baseline_cloud_pct": baseline_raw["meta"]["cloud_cover_pct"],
        "recent_scene": recent_raw["meta"]["scene_id"],
        "recent_date": recent_raw["meta"]["date"],
        "recent_cloud_pct": recent_raw["meta"]["cloud_cover_pct"],
        "n_cells_baseline": int(len(gdf_baseline)),
        "n_cells_recent": int(len(gdf_recent)),
        "mean_lst_baseline_c": round(float(gdf_baseline["LST_C"].mean()), 2),
        "mean_lst_recent_c": round(float(gdf_recent["LST_C"].mean()), 2),
        "warming_delta_c": round(float(warming), 2),
        "urban_mean_lst_c": round(float(urban_mean), 2),
        "green_mean_lst_c": round(float(green_mean), 2),
        "uhi_intensity_c": round(float(uhi_intensity), 2),
        "pearson_r_ndvi_lst": round(float(r_ndvi), 3),
        "pearson_p_ndvi_lst": float(p_ndvi),
        "pearson_r_ndbi_lst": round(float(r_ndbi), 3),
        "pearson_p_ndbi_lst": float(p_ndbi),
        "hotspot_counts": gdf_hotspots["hotspot_class"].value_counts().to_dict(),
    }
    with open(DATA_DIR / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("\n--- Real results summary ---")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print("\nPlotting charts from real data...")

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(gdf_baseline["LST_C"], label=f"{baseline_raw['meta']['date']} (Landsat 8)", fill=True, ax=ax)
    sns.kdeplot(gdf_recent["LST_C"], label=f"{recent_raw['meta']['date']} (Landsat 9)", fill=True, ax=ax)
    ax.set_xlabel("Land Surface Temperature (°C)")
    ax.set_ylabel("Density")
    ax.set_title("Real Tehran LST Distribution — Same Calendar Day, 9 Years Apart")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "lst_histogram.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.regplot(data=valid_pairs, x="NDVI", y="LST_C", scatter_kws={"alpha": 0.35, "s": 10}, line_kws={"color": "red"}, ax=ax)
    ax.set_title(f"NDVI vs. LST (real data) — r = {r_ndvi:.2f}, p < 0.001" if p_ndvi < 0.001 else f"r = {r_ndvi:.2f}")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ndvi_lst_scatter.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1, ax=ax)
    ax.set_title("Real Index Correlation Matrix")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "correlation_heatmap.png", dpi=200)
    plt.close(fig)

    hotspot_colors = {
        "Hot spot - 99% confidence": "#7f0000",
        "Hot spot - 95% confidence": "#d7301f",
        "Hot spot - 90% confidence": "#fc8d59",
        "Not significant": "#e0e0e0",
        "Cold spot - 90% confidence": "#91bfdb",
        "Cold spot - 95% confidence": "#4575b4",
        "Cold spot - 99% confidence": "#313695",
    }
    fig, ax = plt.subplots(figsize=(8, 8))
    for cls, color in hotspot_colors.items():
        subset = gdf_hotspots[gdf_hotspots["hotspot_class"] == cls]
        if len(subset):
            subset.plot(ax=ax, color=color, edgecolor="white", linewidth=0.15, label=cls)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=9)
    ax.set_title("Real Tehran Hot/Cold Clusters — Getis-Ord Gi* (Aug 5, 2024)")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "hotspot_map.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, layer, title in zip(axes, [baseline, recent], [baseline_raw["meta"]["date"], recent_raw["meta"]["date"]]):
        im = ax.imshow(layer["lst_c"], cmap="inferno", vmin=20, vmax=55)
        ax.set_title(f"Real LST — {title}")
        ax.set_axis_off()
    fig.colorbar(im, ax=axes, shrink=0.7, label="°C")
    fig.savefig(FIG_DIR / "lst_raster_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\nAll real outputs saved to {FIG_DIR} and {DATA_DIR}.")


if __name__ == "__main__":
    main()
