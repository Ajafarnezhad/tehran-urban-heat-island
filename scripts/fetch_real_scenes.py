"""دانلود واقعی صحنه‌های Landsat روی تهران از کاتالوگ عمومی Microsoft Planetary Computer.

این اسکریپت جایگزین سریع Earth Engine برای اجرای محلی/بدون احراز هویت گوگل است:
همان الگوریتم‌های src/lst.py و src/indices.py را روی باندهای واقعی دانلودشده اجرا می‌کند
تا نتایج این پروژه با تصاویر واقعی ماهواره‌ای (نه داده شبیه‌سازی‌شده) تولید شوند.

صحنه‌های انتخابی (کمترین ابرناکی، تاریخ تقویمی یکسان برای مقایسه منصفانه):
    baseline : LC08_L2SP_164035_20150805_02_T1  (۵ اوت ۲۰۱۵، ابر ۰.۹٪)
    recent   : LC09_L2SP_164035_20240805_02_T1  (۵ اوت ۲۰۲۴، ابر ۰.۵٪)
"""
import json
from pathlib import Path

import numpy as np
import pystac_client
import planetary_computer
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds

DATA_DIR = Path(__file__).resolve().parent.parent / "outputs" / "data"
RAW_DIR = DATA_DIR / "raw_bands"
RAW_DIR.mkdir(parents=True, exist_ok=True)

TEHRAN_BBOX_WGS84 = [51.15, 35.60, 51.50, 35.80]  # زیرمجموعهٔ فشرده‌تر حول شهر برای دانلود سریع‌تر

SCENES = {
    "baseline": {
        "id": "LC08_L2SP_164035_20150805_02_T1",
        "label": "2015-08-05",
    },
    "recent": {
        "id": "LC09_L2SP_164035_20240805_02_T1",
        "label": "2024-08-05",
    },
}

BANDS = {
    "red": "red",
    "nir08": "nir08",
    "swir16": "swir16",
    "green": "green",
    "lwir11": "lwir11",
    "qa_pixel": "qa_pixel",
}


def get_signed_item(item_id: str):
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace,
    )
    search = catalog.search(collections=["landsat-c2-l2"], ids=[item_id])
    items = list(search.items())
    if not items:
        raise RuntimeError(f"item not found: {item_id}")
    return items[0]


def read_band_window(href: str, bbox_wgs84):
    with rasterio.open(href) as src:
        bounds = transform_bounds("EPSG:4326", src.crs, *bbox_wgs84)
        window = from_bounds(*bounds, transform=src.transform)
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        crs = src.crs
    return data, transform, crs


def fetch_scene(period_key: str):
    scene = SCENES[period_key]
    print(f"[{period_key}] در حال دریافت متادیتای صحنه {scene['id']} ...")
    item = get_signed_item(scene["id"])

    arrays = {}
    transform = None
    crs = None
    for local_name, asset_key in BANDS.items():
        href = item.assets[asset_key].href
        print(f"[{period_key}]   دانلود باند {asset_key} ...")
        data, transform, crs = read_band_window(href, TEHRAN_BBOX_WGS84)
        arrays[local_name] = data

    out_path = RAW_DIR / f"{period_key}.npz"
    np.savez_compressed(
        out_path,
        **arrays,
        transform=np.array(transform)[:6],
        crs=str(crs),
    )
    meta = {
        "scene_id": scene["id"],
        "date": scene["label"],
        "cloud_cover_pct": item.properties.get("eo:cloud_cover"),
        "platform": item.properties.get("platform"),
        "shape": list(next(iter(arrays.values())).shape),
        "crs": str(crs),
    }
    with open(RAW_DIR / f"{period_key}_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    print(f"[{period_key}] ذخیره شد -> {out_path}  (ابر: {meta['cloud_cover_pct']}٪)")
    return meta


def main():
    metas = {}
    for period_key in SCENES:
        metas[period_key] = fetch_scene(period_key)
    print("\nخلاصه صحنه‌های واقعی دانلودشده:")
    print(json.dumps(metas, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
