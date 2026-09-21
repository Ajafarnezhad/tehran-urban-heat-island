"""بارگذاری نقشه پوشش اراضی ESA WorldCover و استخراج ماسک‌های کلاسی برای آمار ناحیه‌ای."""
import ee

WORLDCOVER_CLASSES = {
    10: "جنگل",
    20: "بوته‌زار",
    30: "علفزار",
    40: "زمین کشاورزی",
    50: "مناطق ساخته‌شده",
    60: "پوشش گیاهی/خاک لخت",
    70: "برف/یخ",
    80: "آب",
    90: "تالاب",
    95: "مانگرو",
    100: "خزه/گلسنگ",
}


def load_worldcover(aoi: "ee.Geometry", collection_id: str, year: str = "2021") -> "ee.Image":
    coll = ee.ImageCollection(collection_id)
    image = coll.first() if year == "2021" else coll.filterDate(f"{year}-01-01", f"{year}-12-31").first()
    return image.select("Map").clip(aoi)


def builtup_fraction(worldcover: "ee.Image", region: "ee.Geometry", scale: int = 100) -> float:
    """درصد سطح ساخته‌شده در یک منطقه، برای گزارش سریع در نوت‌بوک."""
    builtup_mask = worldcover.eq(50)
    stats = builtup_mask.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=region, scale=scale, maxPixels=1e9
    )
    return stats.get("Map")
