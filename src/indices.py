"""Spectral index computation from Landsat Collection 2 Level-2 (surface reflectance) imagery."""
import ee


def mask_clouds_l2(image: ee.Image) -> ee.Image:
    """Apply a cloud/shadow mask using the QA_PIXEL band (CFMask bits)."""
    qa = image.select("QA_PIXEL")
    cloud_shadow = 1 << 4
    cloud = 1 << 3
    cirrus = 1 << 2
    mask = (
        qa.bitwiseAnd(cloud_shadow).eq(0)
        .And(qa.bitwiseAnd(cloud).eq(0))
        .And(qa.bitwiseAnd(cirrus).eq(0))
    )
    return image.updateMask(mask)


def apply_scale_factors(image: ee.Image) -> ee.Image:
    """Convert raw L2 optical/thermal band values to reflectance/temperature using official USGS scale factors."""
    optical = image.select("SR_B.").multiply(0.0000275).add(-0.2)
    thermal = image.select("ST_B10").multiply(0.00341802).add(149.0)
    return image.addBands(optical, None, True).addBands(thermal, None, True)


def add_ndvi(image: ee.Image) -> ee.Image:
    ndvi = image.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI")
    return image.addBands(ndvi)


def add_ndbi(image: ee.Image) -> ee.Image:
    """Normalized Difference Built-up Index."""
    ndbi = image.normalizedDifference(["SR_B6", "SR_B5"]).rename("NDBI")
    return image.addBands(ndbi)


def add_ndwi(image: ee.Image) -> ee.Image:
    """McFeeters NDWI, used to exclude water pixels from the analysis."""
    ndwi = image.normalizedDifference(["SR_B3", "SR_B5"]).rename("NDWI")
    return image.addBands(ndwi)


def preprocess_collection(
    collection_id: str,
    aoi: ee.Geometry,
    start: str,
    end: str,
    cloud_cover_max: int = 20,
) -> ee.ImageCollection:
    """Load, cloud-filter, mask and add spectral indices to a Landsat L2 collection."""
    coll = (
        ee.ImageCollection(collection_id)
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt("CLOUD_COVER", cloud_cover_max))
        .map(mask_clouds_l2)
        .map(apply_scale_factors)
        .map(add_ndvi)
        .map(add_ndbi)
        .map(add_ndwi)
    )
    return coll
