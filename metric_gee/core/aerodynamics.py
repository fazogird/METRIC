"""
Module 7: Aerodinamika — shamol, qo'pollik, rah.

F.32: u200 = uw·ln(200/zomw) / ln(zx/zomw) → ERA5 u10 dan
F.33: zom = 0.018·LAI (qisqa ekinlar)
F.34b: zom = exp(a1·NDVI/α + b1) (ko'p turdagi o'simlik)
F.35: zom_mtn = zom·(1 + ((180/π)·s - 5)/20) (tog' tuzatma)
F.36: ω̄ = 1 + 0.1·(elev - elev_station)/1000 (shamol tuzatma)
F.31: u* = k·u200 / ln(200/zom)
F.30: rah = ln(z2/z1) / (u*·k)
F.37: ρair = 1000P / (1.01·(Ts-dT)·R)
"""
import ee
import math
from metric_gee.config.constants import (
    VON_KARMAN, Z_BLEND, Z1, Z2, ZOM_MIN, R_AIR, P_STD, T_STD
)


class Aerodynamics:
    """Shamol va aerodinamik qarshilik hisoblash."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_all(self, era5, dem, lai, ndvi, albedo, ts):
        """Barcha aerodinamik parametrlarni hisoblash."""
        u200 = self.calc_u200(era5, dem)
        zom = self.calc_zom(lai, ndvi, albedo, dem)
        u_star = self.calc_ustar(u200, zom)
        rah = self.calc_rah_neutral(u_star)
        rho_air = self.calc_rho_air(dem, ts)
        
        return {
            'u200': u200,
            'zom': zom,
            'u_star': u_star,
            'rah': rah,
            'rho_air': rho_air,
        }
    
    def calc_u200(self, era5, dem):
        """
        F.32 modified: ERA5 10m shamol → 200m blending height.
        u200 = ws10·ln(200/zomw) / ln(10/zomw)
        
        ERA5 uchun zx = 10m, zomw = 0.03 (o'tloq default).
        """
        ws10 = era5.select('wind_speed_10m')
        zomw = self.settings.zomw
        
        u200 = ws10.multiply(
            math.log(Z_BLEND / zomw) / math.log(10 / zomw)
        ).rename('u200')
        
        # Tog'li shamol tuzatmasi — F.36
        if self.settings.apply_terrain_correction:
            elevation = dem.select('elevation')
            z_mean = elevation.reduceRegion(
                ee.Reducer.mean(), elevation.geometry(), 1000, maxPixels=1e9)
            z_station = ee.Number(z_mean.get('elevation'))
            omega = elevation.subtract(z_station).divide(1000).multiply(0.1).add(1)
            u200 = u200.multiply(omega)
        
        return u200.max(2.0)  # ASCE-EWRI minimal shamol
    
    def calc_zom(self, lai, ndvi, albedo, dem):
        """
        Momentum qo'pollik uzunligi.
        
        F.33: zom = 0.018·LAI (qisqa ekinlar, < 1m)
        F.34b: zom = exp(a1·NDVI/α + b1) (umumiy)
        """
        method = self.settings.zom_method
        
        if method == 'lai':
            # F.33
            zom = lai.multiply(0.018).max(ZOM_MIN)
        elif method == 'ndvi_alpha':
            # F.34b — a1, b1 hududga bog'liq (default qiymatlar)
            a1 = 5.62
            b1 = -5.81
            ratio = ndvi.divide(albedo.max(0.01))
            zom = ratio.multiply(a1).add(b1).exp().max(ZOM_MIN)
        else:
            zom = lai.multiply(0.018).max(ZOM_MIN)
        
        # F.35: Tog'li tuzatma
        if self.settings.apply_terrain_correction:
            slope_deg = dem.select('slope')
            factor = slope_deg.subtract(5).divide(20).add(1).max(1)
            zom = zom.multiply(factor)
        
        return zom.rename('zom')
    
    def calc_ustar(self, u200, zom):
        """F.31: u* = k·u200 / ln(200/zom) — neytral, birinchi iteratsiya."""
        u_star = (u200.multiply(VON_KARMAN)
            .divide(ee.Image.constant(Z_BLEND).divide(zom).log())
            .rename('u_star'))
        return u_star.max(0.01)
    
    def calc_rah_neutral(self, u_star):
        """F.30: rah = ln(z2/z1) / (u*·k) — neytral, birinchi iteratsiya."""
        ln_ratio = math.log(Z2 / Z1)
        rah = u_star.multiply(VON_KARMAN).pow(-1).multiply(ln_ratio).rename('rah')
        return rah
    
    def calc_rho_air(self, dem, ts, dT=None):
        """
        F.37: ρair = 1000P / (1.01·(Ts-dT)·R)
        Birinchi iteratsiyada dT = 0.
        """
        elevation = dem.select('elevation')
        
        # F.5: P = 101.3·((293-0.0065z)/293)^5.26
        P = (elevation.multiply(-0.0065).add(T_STD)
            .divide(T_STD)
            .pow(5.26)
            .multiply(P_STD))
        
        if dT is None:
            dT = ee.Image.constant(0)
        
        # ρair = 1000P / (1.01·(Ts-dT)·R)
        rho = (P.multiply(1000)
            .divide(ts.subtract(dT).multiply(1.01).multiply(R_AIR))
            .rename('rho_air'))
        
        return rho
