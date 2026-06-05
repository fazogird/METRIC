"""
Module 5: Sof radiatsiya (Rn).

F.2: Rn = RS↓ - α·RS↓ + RL↓ - RL↑ - (1-ε₀)·RL↓
F.3: RS↓ = Gsc·cos(θrel)·τsw / d²
F.4: τsw = 0.35 + 0.627·exp[(-0.00146P)/(Kt·cosθhor) - 0.075(W/cosθhor)^0.4]
F.16: RL↑ = ε₀·σ·Ts⁴
F.24: RL↓ = εa·σ·Ta⁴
F.25: εa = 0.85·(-ln(τsw))^0.09

ERA5 dan: P (surface_pressure), ea → W, Ta (temperature_2m)
"""
import ee
from metric_gee.config.constants import G_SC, SIGMA


class NetRadiation:
    """Sof radiatsiya hisoblash."""
    
    def __init__(self, settings):
        self.settings = settings
    
    def calc_rn(self, albedo, ts, e0, era5, dem, solar):
        """
        To'liq Rn hisoblash.
        
        Args:
            albedo: ee.Image — surface albedo
            ts: ee.Image — Ts (K)
            e0: ee.Image — broad-band emissivity
            era5: ee.Image — T2m, ea, sp bands
            dem: ee.Image — elevation
            solar: dict — DOY, d2, cos_theta, sun_elevation
            
        Returns:
            tuple: (rn, tau_sw) — ee.Image
        """
        tau_sw = self._calc_tau_sw(era5, dem, solar)
        rs_down = self._calc_rs_down(tau_sw, solar)
        rl_down = self._calc_rl_down(tau_sw, era5)
        rl_up = self._calc_rl_up(ts, e0)
        
        # F.2: Rn = RS↓ - α·RS↓ + RL↓ - RL↑ - (1-ε₀)·RL↓
        rn = (rs_down
            .subtract(albedo.multiply(rs_down))
            .add(rl_down)
            .subtract(rl_up)
            .subtract(ee.Image.constant(1).subtract(e0).multiply(rl_down))
            .rename('Rn'))
        
        return rn, tau_sw
    
    def _calc_tau_sw(self, era5, dem, solar):
        """
        F.4: Atmosfera quyosh nurini o'tkazuvchanligi.
        τsw = 0.35 + 0.627·exp[(-0.00146P)/(Kt·cosθhor) - 0.075·(W/cosθhor)^0.4]
        """
        # P (kPa) — ERA5 yoki F.5 hisob
        P = era5.select('sp_kpa')
        
        # W (mm) — F.6: W = 0.14·ea·P + 2.1
        ea = era5.select('ea_kpa')
        W = ea.multiply(P).multiply(0.14).add(2.1)
        
        # cos(θhor) — solar zenith
        cos_theta = solar['cos_theta_hor']  # ee.Image yoki ee.Number
        
        Kt = ee.Image.constant(self.settings.kt)
        
        # Beer qonuni asosida
        term1 = P.multiply(-0.00146).divide(Kt.multiply(cos_theta))
        term2 = W.divide(cos_theta).pow(0.4).multiply(-0.075)
        
        tau_sw = term1.add(term2).exp().multiply(0.627).add(0.35).rename('tau_sw')
        
        return tau_sw.clamp(0.1, 0.9)
    
    def _calc_rs_down(self, tau_sw, solar):
        """
        F.3: RS↓ = Gsc·cos(θrel)·τsw / d²
        
        Tekis hudud uchun θrel = θhor (zenithal burchak).
        Tog'li hudud uchun F.7 (slope/aspect) ishlatiladi.
        """
        cos_theta = solar['cos_theta_rel']  # Tog'li: F.7, Tekis: F.8
        d2 = solar['d2']                     # F.9
        
        rs_down = (ee.Image.constant(G_SC)
            .multiply(cos_theta)
            .multiply(tau_sw)
            .divide(d2)
            .rename('RS_down'))
        
        return rs_down.max(0)
    
    def _calc_rl_down(self, tau_sw, era5):
        """
        F.24: RL↓ = εa·σ·Ta⁴
        F.25: εa = 0.85·(-ln(τsw))^0.09
        
        Ta = ERA5 2m temperature (K) — Ts surrogat emas.
        """
        # F.25: Atmosfera emissiviteti
        ea_emiss = tau_sw.log().multiply(-1).pow(0.09).multiply(0.85).rename('ea_emiss')
        
        # Ta — ERA5 dan
        Ta = era5.select('temperature_2m')
        
        # F.24
        rl_down = ea_emiss.multiply(SIGMA).multiply(Ta.pow(4)).rename('RL_down')
        
        return rl_down
    
    def _calc_rl_up(self, ts, e0):
        """
        F.16: RL↑ = ε₀·σ·Ts⁴
        """
        rl_up = e0.multiply(SIGMA).multiply(ts.pow(4)).rename('RL_up')
        return rl_up
