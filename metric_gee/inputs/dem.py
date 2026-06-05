"""
SRTM 30m DEM — balandlik, qiyalik, aspekt.

Terrain korreksiya uchun:
  - TsDEM lapse rate (Mod.3) — F.3 in Allen 2013
  - cos(θrel) qiyalik uchun (Mod.5) — F.7
  - zom_mtn tog'li tuzatma (Mod.7) — F.35
  - Wind speed tuzatma (Mod.7) — F.36
  - Crad slope correction (Mod.11) — F.56
  - Atmospheric pressure (Mod.5) — F.5
"""
import ee
import math


class DEMLoader:
    """SRTM 30m raqamli balandlik modeli."""
    
    COLLECTION_ID = "USGS/SRTMGL1_003"
    
    def __init__(self, settings):
        self.settings = settings
    
    def get_terrain(self, geometry):
        """
        DEM va undan olingan slope/aspect.
        
        Returns:
            ee.Image — 'elevation', 'slope', 'aspect', 'slope_rad', 'aspect_rad' bandlari
        """
        dem = ee.Image(self.COLLECTION_ID).clip(geometry)
        elevation = dem.select('elevation')
        
        # ee.Terrain — slope (degrees), aspect (degrees)
        terrain = ee.Terrain.products(elevation)
        slope_deg = terrain.select('slope')
        aspect_deg = terrain.select('aspect')
        
        # Radian ga o'tkazish (formulalar uchun)
        deg2rad = math.pi / 180
        slope_rad = slope_deg.multiply(deg2rad).rename('slope_rad')
        aspect_rad = aspect_deg.multiply(deg2rad).rename('aspect_rad')
        
        return (elevation
            .addBands(slope_deg)
            .addBands(aspect_deg)
            .addBands(slope_rad)
            .addBands(aspect_rad))
    
    def get_elevation(self, geometry):
        """Faqat balandlik."""
        return ee.Image(self.COLLECTION_ID).select('elevation').clip(geometry)
    
    def get_mean_elevation(self, geometry):
        """AOI o'rtacha balandligi — z_datum uchun (CIMEC)."""
        elevation = self.get_elevation(geometry)
        stats = elevation.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=30,
            maxPixels=1e9
        )
        return ee.Number(stats.get('elevation'))
