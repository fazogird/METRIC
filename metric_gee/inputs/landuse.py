"""
Yer foydalanish maskasi — CIMEC kalibrlash uchun AOI.

Allen 2013: AOI >80% agricultural piksellar bo'lishi kerak.
Cho'l, o'rmon, suv, tik qiyaliklar chiqarilishi kerak.

GEE manbalari:
  - Dynamic World (10m, real-time)
  - ESA WorldCover (10m)
  - USDA CDL (30m, faqat AQSh)
"""
import ee


class LandUseLoader:
    """Qishloq xo'jaligi maskasi generatori."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def get_agri_mask(self, geometry, date_str, source='dynamic_world'):
        """
        Qishloq xo'jaligi binary maskasi — CIMEC AOI uchun.
        
        Args:
            source: 'dynamic_world', 'esa_worldcover', 'cdl'
            
        Returns:
            ee.Image — 1 = agricultural, 0 = boshqa
        """
        if source == 'dynamic_world':
            return self._mask_dynamic_world(geometry, date_str)
        elif source == 'esa_worldcover':
            return self._mask_esa(geometry)
        elif source == 'cdl':
            return self._mask_cdl(geometry, date_str)
        else:
            raise ValueError(f"Noma'lum manba: {source}")
    
    def _mask_dynamic_world(self, geometry, date_str):
        """
        Dynamic World — 'crops' klassi.
        Global qamrov, 10m, Sentinel-2 asosida.
        """
        date = ee.Date(date_str)
        
        dw = (ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
            .filterBounds(geometry)
            .filterDate(date.advance(-30, 'day'), date.advance(30, 'day'))
            .select('label')
            .mode())  # Eng ko'p uchraydigan klass
        
        # Dynamic World labels: 4 = crops
        crops_mask = dw.eq(4).rename('agri_mask')
        
        return crops_mask.clip(geometry)
    
    def _mask_esa(self, geometry):
        """ESA WorldCover 2021 — 10m global."""
        esa = ee.Image("ESA/WorldCover/v200/2021")
        
        # 40 = Cropland
        crops_mask = esa.eq(40).rename('agri_mask')
        
        return crops_mask.clip(geometry)
    
    def _mask_cdl(self, geometry, date_str):
        """USDA Cropland Data Layer — 30m, faqat AQSh."""
        year = ee.Date(date_str).get('year')
        
        cdl = (ee.ImageCollection('USDA/NASS/CDL')
            .filter(ee.Filter.calendarRange(year, year, 'year'))
            .first()
            .select('cropland'))
        
        # CDL: 1-60 = row crops, 61-77 = boshqa ekinlar
        crops_mask = cdl.lte(77).And(cdl.gte(1)).rename('agri_mask')
        
        return crops_mask.clip(geometry)
    
    def get_water_mask(self, geometry):
        """Suv maskasi — G/Rn = 0.5 uchun."""
        water = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
        return water.select('occurrence').gt(50).rename('water_mask').clip(geometry)
