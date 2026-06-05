"""
METRIC-GEE loyiha sozlamalari.
Har bir run uchun o'zgartiriladi.
"""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Settings:
    """Pipeline konfiguratsiyasi."""
    
    # --- Study Area ---
    study_area_name: str = "Idaho_2024"
    
    # --- Sana oralig'i ---
    start_date: str = "2024-06-01"
    end_date: str = "2024-09-30"
    
    # --- Landsat ---
    cloud_cover_max: int = 20          # % bulut filtri
    landsat_collection: str = "LANDSAT/LC08/C02/T1_L2"  # yoki LC09
    
    # --- ERA5 ---
    era5_collection: str = "ECMWF/ERA5_LAND/HOURLY"
    
    # --- DEM ---
    dem_collection: str = "USGS/SRTMGL1_003"
    
    # --- Land Use ---
    landuse_collection: str = "GOOGLE/DYNAMICWORLD/V1"
    
    # --- Albedo ---
    albedo_method: str = "olmedo"      # 'olmedo' yoki 'liang'
    
    # --- Soil Heat Flux ---
    g_method: str = "bastiaanssen"     # 'bastiaanssen' (F.26) yoki 'tasumi' (F.27)
    
    # --- Roughness ---
    zom_method: str = "lai"            # 'lai' (F.33) yoki 'ndvi_alpha' (F.34b)
    
    # --- CIMEC Anchor ---
    use_auto_calibration: bool = True  # False → manual pixel tanlash
    lapse_rate: float = 6.5            # K/km

    # ESA WorldCover v100 LULC kodlar:
    #   10-Trees  20-Shrub  30-Grass  40-Crop  50-Built  60-Bare  80-Water
    # Cold piksel (sug'oriladigan ekin):
    cold_lulc_codes: list = None       # None → [40] (Cropland)
    # Hot piksel (yalang'och tuproq / shahar):
    hot_lulc_codes:  list = None       # None → [60, 50] (Bare + Built-up)
    
    # --- Terrain ---
    apply_terrain_correction: bool = True  # F.7 slope/aspect, F.35-36
    
    # --- Meteorologiya ---
    kt: float = 1.0                    # Turbidity coefficient (F.4)
    zomw: float = 0.03                 # Meteostansiya roughness (o'tloq)
    timezone_lon_deg: float = 0.0      # Vaqt zonasi: Idaho=-105, Uz=+69
    
    # --- FAO56 Bare Soil ---
    etrf_bare: float = 0.05            # Default bare soil ETrF
    
    # --- Export ---
    export_scale: int = 30             # m
    export_crs: str = "EPSG:4326"
    export_folder: str = "METRIC_outputs"
    
    # --- Landsat tile filter (validatsiya uchun bir tile) ---
    wrs_path: Optional[int] = None   # WRS-2 Path (masalan Idaho: 39)
    wrs_row:  Optional[int] = None   # WRS-2 Row  (masalan Idaho: 29)

    # --- Validation ---
    openet_collection: Optional[str] = "OpenET/ENSEMBLE/CONUS/GRIDMET/MONTHLY/v2_0"
    
    def get_landsat_collections(self) -> list:
        """L8 va L9 ikkala collectionni qaytaradi."""
        return [
            "LANDSAT/LC08/C02/T1_L2",
            "LANDSAT/LC09/C02/T1_L2"
        ]
