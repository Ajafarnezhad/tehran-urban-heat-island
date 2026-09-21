"""راه‌اندازی اتصال به Google Earth Engine و بارگذاری پیکربندی پروژه."""
import yaml
import ee
from pathlib import Path

CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "config.yaml"


def load_config(path: Path = CONFIG_PATH) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def init_earth_engine(project_id: str | None = None) -> None:
    """احراز هویت و مقداردهی اولیه Earth Engine.

    اولین اجرا مرورگر را برای ورود با حساب گوگل باز می‌کند؛ توکن حاصل
    به صورت محلی کش می‌شود و اجراهای بعدی نیازی به ورود دوباره ندارند.
    """
    cfg = load_config()
    project_id = project_id or cfg["project"]["ee_project_id"]
    try:
        ee.Initialize(project=project_id)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project_id)


def get_aoi(cfg: dict) -> ee.Geometry:
    min_lon, min_lat, max_lon, max_lat = cfg["aoi"]["bbox"]
    return ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
