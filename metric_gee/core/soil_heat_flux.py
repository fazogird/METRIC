"""
Module 6: Tuproq issiqlik oqimi (G).

F.26 (Bastiaanssen 2000):
  G/Rn = (Ts-273.15)(0.0038 + 0.0074α)(1 - 0.98·NDVI⁴)
  
F.27a (Tasumi, LAI ≥ 0.5):
  G/Rn = 0.05 + 0.18·exp(-0.521·LAI)

F.27b (Tasumi, LAI < 0.5):
  G/Rn = 1.80·(Ts-273.15)/Rn + 0.084

Suv: G/Rn = 0.5
Qor: G/Rn = 0.5 (Ts < 277K, NDVI < 0)
"""
import ee


class SoilHeatFlux:
    """Tuproq issiqlik oqimi."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_g(self, rn, ts, ndvi, albedo, lai=None):
        """
        G hisoblash — tanlangan usul bilan.
        
        Returns:
            ee.Image — 'G' (W/m²)
        """
        method = self.settings.g_method
        
        if method == 'bastiaanssen':
            g_ratio = self._bastiaanssen(ts, ndvi, albedo)
        elif method == 'tasumi':
            g_ratio = self._tasumi(ts, rn, lai)
        else:
            raise ValueError(f"Noma'lum G usuli: {method}")
        
        # Suv va qor uchun maxsus qiymatlar
        g_ratio = g_ratio.where(ndvi.lt(0).And(albedo.lt(0.1)), 0.5)   # Suv
        g_ratio = g_ratio.where(ts.lt(277).And(ndvi.lt(0)), 0.5)        # Qor
        
        g = g_ratio.multiply(rn).rename('G')
        return g
    
    def _bastiaanssen(self, ts, ndvi, albedo):
        """
        F.26: G/Rn = (Ts-273.15)(0.0038 + 0.0074α)(1 - 0.98·NDVI⁴)
        
        Universal formula — global model uchun default.
        """
        ts_c = ts.subtract(273.15)
        
        g_ratio = (ts_c
            .multiply(albedo.multiply(0.0074).add(0.0038))
            .multiply(ee.Image.constant(1).subtract(ndvi.pow(4).multiply(0.98)))
            .rename('G_ratio'))
        
        return g_ratio.clamp(0, 0.65)
    
    def _tasumi(self, ts, rn, lai):
        """
        F.27a: G/Rn = 0.05 + 0.18·exp(-0.521·LAI)  (LAI ≥ 0.5)
        F.27b: G/Rn = 1.80·(Ts-273.15)/Rn + 0.084   (LAI < 0.5)
        
        Idaho-specific, Wright (1982) ma'lumotlari asosida.
        """
        ts_c = ts.subtract(273.15)
        
        # F.27a: o'simlikli
        g_veg = lai.multiply(-0.521).exp().multiply(0.18).add(0.05)
        
        # F.27b: yalang'och tuproq
        g_bare = ts_c.multiply(1.80).divide(rn).add(0.084)
        
        g_ratio = g_veg.where(lai.lt(0.5), g_bare).rename('G_ratio')
        
        return g_ratio.clamp(0, 0.65)
