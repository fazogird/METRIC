"""QA_PIXEL bitwise cloud/shadow mask — Landsat 8/9 C2."""
import ee


def apply_cloud_mask(image):
    """QA_PIXEL dan cloud + shadow mask."""
    qa = image.select('QA_PIXEL')
    cloud = qa.bitwiseAnd(1 << 4).eq(0)
    shadow = qa.bitwiseAnd(1 << 3).eq(0)
    return image.updateMask(cloud.And(shadow))
