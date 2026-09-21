"""Export final layers to Google Drive for archiving and use in desktop GIS."""
import ee


def export_image_to_drive(image: ee.Image, description: str, region: ee.Geometry, cfg: dict):
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=description,
        folder=cfg["export"]["drive_folder"],
        region=region,
        scale=cfg["export"]["scale_m"],
        crs=cfg["export"]["crs"],
        maxPixels=1e10,
        fileFormat="GeoTIFF",
    )
    task.start()
    return task


def export_table_to_drive(feature_collection: ee.FeatureCollection, description: str, cfg: dict):
    task = ee.batch.Export.table.toDrive(
        collection=feature_collection,
        description=description,
        folder=cfg["export"]["drive_folder"],
        fileFormat="CSV",
    )
    task.start()
    return task
