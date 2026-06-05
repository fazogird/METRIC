"""
Module 3: Yer yuzasi temperaturasi (Ts) va TsDEM.

F.20-21 SKIP — L2 ST_B10 allaqachon USGS TICOS orqali tuzatilgan.
  Rp, τNB, Rsky — kerak emas.

TsDEM — Allen 2013 (CIMEC uchun):
  TsDEM = Ts - lapse_rate × (z - z_datum) / 1000
  
  Lapse rate: 6.5 K/km default (troposfera standarti).
  z_datum: AOI o'rtacha balandligi.
"""
import ee
from metric_gee.config.constants import LAPSE_RATE_DEFAULT


class SurfaceTemperature:
    """Ts va TsDEM hisoblash."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_ts(self, image):
        """
        L2 Surface Temperature — allaqachon Kelvin da.
        Scale factor pipeline.py da qo'llanilgan (LandsatLoader).
        """
        return image.select('B_THERMAL').rename('Ts')
    
    def calc_ts_dem(self, ts, dem, z_datum=None):
        """
        TsDEM — balandlik effektini olib tashlash.
        
        TsDEM = Ts - lapse_rate × (z - z_datum) / 1000
        
        Bu qilinadi chunki baland joylarda Ts tabiiyan past —
        lekin bu bug'lanish tufayli emas, atmosfera sovuqligi tufayli.
        TsDEM faqat bug'lanish ta'sirini ko'rsatadi.
        
        Args:
            ts: ee.Image — Ts (K)
            dem: ee.Image — 'elevation' band (m)
            z_datum: ee.Number — AOI o'rtacha balandlik (m)
            
        Returns:
            ee.Image — 'TsDEM' (K)
        """
        elevation = dem.select('elevation')
        lapse_rate = self.settings.lapse_rate  # K/km
        
        # z_datum = AOI o'rtacha balandligi (agar berilmagan bo'lsa)
        if z_datum is None:
            z_datum = ee.Number(0)  # dengiz sathi → pipeline.py da to'g'ri qiymat beriladi
        
        # Balandlik farqi (m) → km ga o'tkazish
        delta_z = elevation.subtract(z_datum).divide(1000)
        
        # TsDEM = Ts - lapse_rate × Δz
        ts_dem = ts.subtract(delta_z.multiply(lapse_rate)).rename('TsDEM')
        
        return ts_dem
