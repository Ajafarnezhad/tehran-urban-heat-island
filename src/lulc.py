"""Load the ESA WorldCover land-cover map and extract class masks for zonal statistics."""
import ee

WORLDCOVER_CLASSES = {
    10: "Tree cover",
    20: "Shrubland",
    30: "Grassland",
    40: "Cropland",
    50: "Built-up",
    60: "Bare / sparse vegetation",
    70: "Snow and ice",
    80: "Permanent water bodies",
    90: "Herbaceous wetland",
    95: "Mangroves",
    100: "Moss and lichen",
}


def load_worldcover(aoi: "ee.Geometry", collection_id: str, year: str = "2021") -> "ee.Image":
    coll = ee.ImageCollection(collection_id)
    image = coll.first() if year == "2021" else coll.filterDate(f"{year}-01-01", f"{year}-12-31").first()
    return image.select("Map").clip(aoi)


def builtup_fraction(worldcover: "ee.Image", region: "ee.Geometry", scale: int = 100) -> float:
    """Fraction of built-up surface within a region, for a quick notebook-level summary."""
    builtup_mask = worldcover.eq(50)
    stats = builtup_mask.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=region, scale=scale, maxPixels=1e9
    )
    return stats.get("Map")
