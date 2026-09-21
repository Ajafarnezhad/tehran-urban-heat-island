"""استخراج دمای سطح زمین (LST) با روش تابش‌سنجی تک‌کاناله (Mono-Window) بر پایه NDVI Emissivity."""
import ee


def add_emissivity(image: ee.Image) -> ee.Image:
    """برآورد گسیل‌مندی سطح (LSE) از روی NDVI با روش آستانه NDVI (NDVI Threshold Method).

    - آب (NDVI < 0): گسیل‌مندی ثابت 0.991
    - خاک لخت (NDVI < 0.2): گسیل‌مندی ثابت 0.966
    - پوشش گیاهی کامل (NDVI > 0.5): گسیل‌مندی ثابت 0.973
    - پیکسل مخلوط: با استفاده از سهم پوشش گیاهی (Pv) و ضریب اثر هندسی سطح (0.005)
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
    """محاسبه LST بر حسب سلسیوس از دمای درخشندگی باند حرارتی (ST_B10، از پیش کالیبره‌شده به کلوین)
    با تصحیح گسیل‌مندی طبق معادلهٔ پلانک ساده‌شده (Artis & Carnahan, 1982؛ روش پرکاربرد USGS L2 LST).
    """
    image = add_emissivity(image)
    lst_kelvin = image.select("ST_B10")
    emissivity = image.select("EM")

    wavelength = 10.895  # میکرومتر، میانگین باند حرارتی TIRS
    rho = 14388.0  # h*c/k_B بر حسب میکرومتر-کلوین (h*c/sigma_B)

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
    """میانگین‌گیری زمانی از یک کالکشن پیش‌پردازش‌شده برای تولید ترکیب LST/NDVI/NDBI فصلی."""
    with_lst = collection.map(add_lst_celsius)
    composite = with_lst.select(["LST_C", "NDVI", "NDBI", "NDWI"]).mean()
    return composite
