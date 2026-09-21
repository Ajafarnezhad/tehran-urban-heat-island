"""نقطه ورود CLI: اجرای کامل پایپ‌لاین تحلیل UHI تهران از ابتدا تا انتها.

اجرا:
    python scripts/run_pipeline.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import ee
import pandas as pd

from src.gee_setup import init_earth_engine, load_config, get_aoi
from src.indices import preprocess_collection
from src.lst import build_lst_composite
from src.lulc import load_worldcover
from src.uhi_analysis import (
    sample_grid_to_dataframe,
    to_geodataframe,
    correlation_report,
    getis_ord_hotspots,
    uhi_intensity,
)
from src.visualization import (
    plot_lst_histogram,
    plot_ndvi_lst_scatter,
    plot_hotspot_map,
    plot_correlation_heatmap,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "outputs"


def build_period_composite(cfg, aoi, period_key, collection_id):
    period = cfg["time_periods"][period_key]
    coll = preprocess_collection(
        collection_id,
        aoi,
        period["start"],
        period["end"],
        cfg["landsat"]["cloud_cover_max"],
    )
    return build_lst_composite(coll)


def main():
    cfg = load_config()
    init_earth_engine(cfg["project"]["ee_project_id"])
    aoi = get_aoi(cfg)

    print("در حال ساخت ترکیب دمایی دوره پایه (2015)...")
    lst_baseline = build_period_composite(cfg, aoi, "baseline", cfg["landsat"]["baseline_collection"])

    print("در حال ساخت ترکیب دمایی دوره اخیر (2024)...")
    lst_recent = build_period_composite(cfg, aoi, "recent", cfg["landsat"]["recent_collection"])

    print("در حال بارگذاری نقشه پوشش اراضی...")
    worldcover = load_worldcover(aoi, cfg["lulc"]["collection"], cfg["lulc"]["year"])

    cell_size = cfg["analysis"]["grid_cell_size_m"]

    print("در حال نمونه‌برداری شبکه‌ای برای دوره پایه...")
    df_baseline = sample_grid_to_dataframe(lst_baseline, aoi, cell_size)
    gdf_baseline = to_geodataframe(df_baseline)

    print("در حال نمونه‌برداری شبکه‌ای برای دوره اخیر...")
    df_recent = sample_grid_to_dataframe(lst_recent, aoi, cell_size)
    gdf_recent = to_geodataframe(df_recent)

    gdf_recent.to_file(OUT_DIR / "data" / "grid_recent.geojson", driver="GeoJSON")
    gdf_baseline.to_file(OUT_DIR / "data" / "grid_baseline.geojson", driver="GeoJSON")

    print("در حال محاسبه ماتریس همبستگی...")
    corr = correlation_report(gdf_recent)
    corr.to_csv(OUT_DIR / "data" / "correlation_matrix.csv", encoding="utf-8-sig")

    print("در حال اجرای تحلیل هات‌اسپات Getis-Ord Gi*...")
    gdf_hotspots = getis_ord_hotspots(gdf_recent, value_col="LST_C")
    gdf_hotspots.to_file(OUT_DIR / "data" / "hotspots.geojson", driver="GeoJSON")

    print("در حال رسم نمودارها...")
    plot_lst_histogram(gdf_baseline, gdf_recent, str(OUT_DIR / "figures" / "lst_histogram.png"))
    plot_ndvi_lst_scatter(gdf_recent, str(OUT_DIR / "figures" / "ndvi_lst_scatter.png"))
    plot_hotspot_map(gdf_hotspots, str(OUT_DIR / "figures" / "hotspot_map.png"))
    plot_correlation_heatmap(corr, str(OUT_DIR / "figures" / "correlation_heatmap.png"))

    summary = {
        "baseline_mean_lst_c": round(gdf_baseline["LST_C"].mean(), 2),
        "recent_mean_lst_c": round(gdf_recent["LST_C"].mean(), 2),
        "warming_delta_c": round(gdf_recent["LST_C"].mean() - gdf_baseline["LST_C"].mean(), 2),
        "correlation_ndvi_lst": round(corr.loc["NDVI", "LST_C"], 3),
    }
    pd.Series(summary).to_csv(OUT_DIR / "data" / "summary.csv", encoding="utf-8-sig")

    print("\n--- خلاصه نتایج ---")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nتمام خروجی‌ها در {OUT_DIR} ذخیره شدند.")


if __name__ == "__main__":
    main()
