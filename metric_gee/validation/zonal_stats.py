"""Zonal statistics — field/district darajasida ET statistika."""
import ee


def calc_field_stats(et_image, fields_fc, scale=30):
    """Har bir field uchun o'rtacha ET."""
    return et_image.reduceRegions(
        collection=fields_fc,
        reducer=ee.Reducer.mean().combine(
            ee.Reducer.stdDev(), sharedInputs=True),
        scale=scale)


def calc_regional_stats(et_image, geometry, scale=30):
    """Hudud bo'yicha umumiy statistika."""
    return et_image.reduceRegion(
        reducer=ee.Reducer.mean().combine(
            ee.Reducer.stdDev(), sharedInputs=True).combine(
            ee.Reducer.minMax(), sharedInputs=True),
        geometry=geometry,
        scale=scale,
        maxPixels=1e9)
