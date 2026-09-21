"""Land Surface Temperature (LST) retrieval using the mono-window method with NDVI-based emissivity."""
import ee


def add_emissivity(image: ee.Image) -> ee.Image:
    """Estimate land surface emissivity (LSE) from NDVI using the NDVI Threshold Method.

    - Water (NDVI < 0): constant emissivity 0.991
    - Bare soil (NDVI < 0.2): constant emissivity 0.966
    - Full vegetation cover (NDVI > 0.5): constant emissivity 0.973
    - Mixed pixel: interpolated from vegetation fraction (Pv) plus a surface-roughness
      correction term (0.005)
    """
    ndvi = image.select("NDVI")
    pv = ndvi.subtract(0.2).divide(0.3).pow(2).clamp(0, 1).rename("PV")

    emissivity_mixed = pv.multiply(0.973 - 0.966).add(0.966).add(0.005).rename("EM")

    emissivity = (
        ee.Image(0.991)
        .where(ndvi.gte(-1).And(ndvi.lt(0)), 0.991)
        .where(ndvi.gte(0).And(ndvi.lt(0.2)), 0.966)
        .where(ndvi.gte(0.2).And(ndvi.lte(0.5)), emissivity_mixed)
        .where(ndvi.gt(0.5), 0.973)
        .rename("EM")
    )
    return image.addBands(emissivity).addBands(pv)


def add_lst_celsius(image: ee.Image) -> ee.Image:
    """Compute LST in Celsius from the thermal band's brightness temperature (ST_B10, already
    calibrated to Kelvin) with an emissivity correction using the simplified Planck equation
    (Artis & Carnahan, 1982; a common approach for USGS L2 LST products).
    """
    image = add_emissivity(image)
    lst_kelvin = image.select("ST_B10")
    emissivity = image.select("EM")

    wavelength = 10.895  # micrometers, mean effective wavelength of the TIRS thermal band
    rho = 14388.0  # h*c/k_B in micrometer-Kelvin

    lst_corrected = lst_kelvin.expression(
        "LST / (1 + (WL * LST / RHO) * log(EM))",
        {
            "LST": lst_kelvin,
            "WL": wavelength,
            "RHO": rho,
            "EM": emissivity,
        },
    ).rename("LST_K")

    lst_celsius = lst_corrected.subtract(273.15).rename("LST_C")
    return image.addBands(lst_corrected).addBands(lst_celsius)


def build_lst_composite(collection) -> "ee.Image":
    """Temporally average a preprocessed collection to produce a seasonal LST/NDVI/NDBI composite."""
    with_lst = collection.map(add_lst_celsius)
    composite = with_lst.select(["LST_C", "NDVI", "NDBI", "NDWI"]).mean()
    return composite
