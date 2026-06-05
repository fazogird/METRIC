"""
Module 2: Yer yuzasi albedosi.

METRIC 2007 (F.10-15): TOA → atmospheric correction → surface albedo
Biz: L2 SR allaqachon surface reflectance → to'g'ridan vaznli yig'indi.

F.10-14 SKIP — L2 LaSRC tuzatilgan.
F.15 MODIFIED — Olmedo yoki Liang vaznlari.

Olmedo et al. — Landsat 8 OLI surface reflectance uchun:
  α = 0.4739·B2 - 0.4372·B3 + 0.1652·B4 + 0.2831·B5 + 0.1072·B6 + 0.1029·B7 + 0.0366

Liang (2001) — Landsat uchun alternativ:
  α = 0.356·B2 + 0.130·B4 + 0.373·B5 + 0.085·B6 + 0.072·B7 - 0.0018

Izoh: Albedo [0, 1] oralig'ida bo'lishi kerak.
  Qor: ~0.8, Suv: ~0.06, O'simlik: 0.15-0.25, Cho'l: 0.30-0.45
"""
import ee


class AlbedoCalculator:
    """Surface albedo hisoblash — L8/9 C2 L2 SR uchun."""
    
    # Olmedo vaznlari — Landsat 8/9 OLI
    OLMEDO_WEIGHTS = {
        'B_BLUE':  0.4739,
        'B_GREEN': -0.4372,
        'B_RED':   0.1652,
        'B_NIR':   0.2831,
        'B_SWIR1': 0.1072,
        'B_SWIR2': 0.1029,
    }
    OLMEDO_OFFSET = 0.0366
    
    # Liang (2001) vaznlari
    LIANG_WEIGHTS = {
        'B_BLUE':  0.356,
        'B_RED':   0.130,
        'B_NIR':   0.373,
        'B_SWIR1': 0.085,
        'B_SWIR2': 0.072,
    }
    LIANG_OFFSET = -0.0018
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_albedo(self, image, method=None):
        """
        Surface albedo hisoblash.
        
        Args:
            image: ee.Image — scaled SR bands (B_BLUE, B_GREEN, ...)
            method: str — 'olmedo' yoki 'liang' (default: settings dan)
            
        Returns:
            ee.Image — 'albedo' band, [0, 1] clamped
        """
        method = method or self.settings.albedo_method
        
        if method == 'olmedo':
            albedo = self._olmedo(image)
        elif method == 'liang':
            albedo = self._liang(image)
        else:
            raise ValueError(f"Noma'lum albedo usuli: {method}")
        
        return albedo.clamp(0, 1).rename('albedo')
    
    def _olmedo(self, image):
        """
        Olmedo et al. — L8/9 OLI SR uchun.
        α = 0.4739·B2 - 0.4372·B3 + 0.1652·B4 + 0.2831·B5 + 0.1072·B6 + 0.1029·B7 + 0.0366
        """
        albedo = ee.Image.constant(self.OLMEDO_OFFSET)
        
        for band, weight in self.OLMEDO_WEIGHTS.items():
            albedo = albedo.add(image.select(band).multiply(weight))
        
        return albedo
    
    def _liang(self, image):
        """
        Liang (2001) — L8/9 uchun.
        α = 0.356·B2 + 0.130·B4 + 0.373·B5 + 0.085·B6 + 0.072·B7 - 0.0018
        """
        albedo = ee.Image.constant(self.LIANG_OFFSET)
        
        for band, weight in self.LIANG_WEIGHTS.items():
            albedo = albedo.add(image.select(band).multiply(weight))
        
        return albedo
