"""
Module 1: O'simlik indekslari — NDVI, SAVI, LAI.

Formulalar (Allen 2007):
  F.23: NDVI = (NIR - RED) / (NIR + RED)
  F.19: SAVI = (1+L)(NIR - RED) / (L + NIR + RED)
  F.18: LAI  = -ln((0.69 - SAVI_ID) / 0.59) / 0.91

L8/9 band mapping:
  RED = B4 (SR_B4) → B_RED
  NIR = B5 (SR_B5) → B_NIR

Izoh: L2 surface reflectance ustida hisob — TOA emas.
"""
import ee


class VegetationIndices:
    """NDVI, SAVI, LAI hisoblash."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_all(self, image):
        """Barcha indekslarni hisoblash."""
        ndvi = self.calc_ndvi(image)
        savi = self.calc_savi(image)
        lai = self.calc_lai(savi)
        return ndvi, savi, lai
    
    def calc_ndvi(self, image):
        """
        F.23: NDVI = (B5 - B4) / (B5 + B4)
        
        L8/9 SR surface reflectance ustida.
        Diapazon: [-1, 1], odatda o'simlik: 0.2–0.9
        """
        nir = image.select('B_NIR')
        red = image.select('B_RED')
        
        ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI')
        return ndvi.clamp(-1, 1)
    
    def calc_savi(self, image, L=0.1):
        """
        F.19: SAVI = (1+L)(NIR - RED) / (L + NIR + RED)
        
        L = 0.1 (Idaho, Tasumi 2003) — tuproq fon ta'sirini kamaytiradi.
        L = 0.5 — umumiy qiymat (Huete 1988).
        
        METRIC da L=0.1 ishlatiladi → SAVI_ID deb ataladi.
        """
        nir = image.select('B_NIR')
        red = image.select('B_RED')
        
        savi = (nir.subtract(red)
            .multiply(1 + L)
            .divide(nir.add(red).add(L))
            .rename('SAVI'))
        
        return savi
    
    def calc_lai(self, savi):
        """
        F.18: LAI = -ln((0.69 - SAVI) / 0.59) / 0.91
        
        Chegaralar:
          LAI = 0   agar SAVI < 0.1
          LAI = 6.0 agar SAVI > 0.687
          
        LAI o'simlik massasini ifodalaydi (m²/m²).
        """
        # Asosiy hisob
        lai = (ee.Image.constant(0.69)
            .subtract(savi)
            .divide(0.59)
            .log()
            .multiply(-1)
            .divide(0.91)
            .rename('LAI'))
        
        # Chegaralar qo'llash
        lai = lai.where(savi.lt(0.1), 0)
        lai = lai.where(savi.gt(0.687), 6.0)
        lai = lai.clamp(0, 6.0)
        
        return lai
