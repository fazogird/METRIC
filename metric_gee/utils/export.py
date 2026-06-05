"""GeoTIFF, Asset, Drive export funksiyalari."""
import ee


def to_drive(image, description, folder, geometry, scale=30):
    """Google Drive ga export."""
    task = ee.batch.Export.image.toDrive(
        image=image, description=description,
        folder=folder, region=geometry,
        scale=scale, maxPixels=1e13,
        fileFormat='GeoTIFF')
    task.start()
    return task


def to_asset(image, asset_id, geometry, scale=30):
    """GEE Asset ga export."""
    task = ee.batch.Export.image.toAsset(
        image=image, description=asset_id.split('/')[-1],
        assetId=asset_id, region=geometry,
        scale=scale, maxPixels=1e13)
    task.start()
    return task
