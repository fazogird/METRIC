"""
Landsat 8/9 Collection 2 Level 2 yuklagich.

L2 = USGS tomonidan atmosferik korreksiya qilingan:
  - SR bands: LaSRC algoritm (surface reflectance)
  - ST band: TICOS algoritm (surface temperature)

Band mapping L5/7 → L8/9:
  L5/7: B1(blue) B2(green) B3(red) B4(nir) B5(swir1) B7(swir2) B6(thermal)
  L8/9: B2(blue) B3(green) B4(red) B5(nir) B6(swir1) B7(swir2) B10(thermal)
"""
import ee
from metric_gee.config.constants import L8_SR_SCALE, L8_SR_OFFSET, L8_ST_SCALE, L8_ST_OFFSET


class LandsatLoader:
    """Landsat 8/9 C2 L2 yuklash va tayyorlash."""
    
    # L8/9 OLI band nomlari → umumiy nomlar
    BAND_MAP = {
        'SR_B2': 'B_BLUE',
        'SR_B3': 'B_GREEN',
        'SR_B4': 'B_RED',
        'SR_B5': 'B_NIR',
        'SR_B6': 'B_SWIR1',
        'SR_B7': 'B_SWIR2',
        'ST_B10': 'B_THERMAL',
        'QA_PIXEL': 'QA_PIXEL',
    }
    
    def __init__(self, settings):
        self.settings = settings
    
    def get_image(self, geometry, date):
        """
        Berilgan sana va hududga eng yaqin Landsat tasvirni topadi.
        
        Args:
            geometry: ee.Geometry — study area
            date: str — 'YYYY-MM-DD'
            
        Returns:
            ee.Image — scaled, renamed, masked
        """
        collections = self.settings.get_landsat_collections()
        
        # L8 va L9 ni birlashtirish
        merged = ee.ImageCollection([])
        for coll_id in collections:
            coll = (ee.ImageCollection(coll_id)
                .filterBounds(geometry)
                .filterDate(date, ee.Date(date).advance(16, 'day'))
                .filter(ee.Filter.lt('CLOUD_COVER', self.settings.cloud_cover_max)))
            merged = merged.merge(coll)
        
        # Eng kam bulutli tasvirni tanlash
        image = merged.sort('CLOUD_COVER').first()
        
        # Qayta ishlash
        image = self._apply_scale_factors(image)
        image = self._rename_bands(image)
        image = self._apply_cloud_mask(image)
        image = self._add_metadata(image)
        
        return image
    
    def get_collection(self, geometry, start_date, end_date):
        """Sana oralig'idagi barcha tasvirlar kolleksiyasi."""
        collections = self.settings.get_landsat_collections()
        
        merged = ee.ImageCollection([])
        for coll_id in collections:
            coll = (ee.ImageCollection(coll_id)
                .filterBounds(geometry)
                .filterDate(start_date, end_date)
                .filter(ee.Filter.lt('CLOUD_COVER', self.settings.cloud_cover_max)))
            merged = merged.merge(coll)
        
        return merged.map(lambda img: self._process_image(img))
    
    def _process_image(self, image):
        """Bitta tasvirni to'liq qayta ishlash."""
        image = self._apply_scale_factors(image)
        image = self._rename_bands(image)
        image = self._apply_cloud_mask(image)
        image = self._add_metadata(image)
        return image
    
    def _apply_scale_factors(self, image):
        """
        L2 scale factors qo'llash.
        SR: reflectance = DN * 0.0000275 - 0.2  → [0, 1] oraliq
        ST: kelvin = DN * 0.00341802 + 149.0
        """
        sr_bands = image.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'])
        sr_scaled = sr_bands.multiply(L8_SR_SCALE).add(L8_SR_OFFSET).clamp(0, 1)
        
        st_band = image.select(['ST_B10'])
        st_scaled = st_band.multiply(L8_ST_SCALE).add(L8_ST_OFFSET)
        
        qa_band = image.select(['QA_PIXEL'])
        
        return sr_scaled.addBands(st_scaled).addBands(qa_band).copyProperties(image, image.propertyNames())
    
    def _rename_bands(self, image):
        """Band nomlarini umumiy nomlarga o'zgartirish."""
        old_names = list(self.BAND_MAP.keys())
        new_names = list(self.BAND_MAP.values())
        return image.select(old_names, new_names)
    
    def _apply_cloud_mask(self, image):
        """QA_PIXEL bitwise cloud/shadow mask."""
        qa = image.select('QA_PIXEL')
        
        # Bit 3 = cloud shadow, Bit 4 = cloud
        cloud_shadow = qa.bitwiseAnd(1 << 3).eq(0)
        cloud = qa.bitwiseAnd(1 << 4).eq(0)
        
        mask = cloud_shadow.And(cloud)
        return image.updateMask(mask)
    
    def _add_metadata(self, image):
        """Kerakli metadata propertylarni qo'shish."""
        # Sun elevation → zenith
        sun_elev = image.get('SUN_ELEVATION')
        
        return (image
            .set('sun_elevation', sun_elev)
            .set('acquisition_time', image.get('SCENE_CENTER_TIME'))
            .set('spacecraft_id', image.get('SPACECRAFT_ID'))
            .set('doy', ee.Date(image.get('system:time_start')).getRelative('day', 'year').add(1)))
