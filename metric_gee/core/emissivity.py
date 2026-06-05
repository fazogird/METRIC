"""
Module 4: Emissivitet — ε₀ (broad-band) va εNB (narrow-band).

F.17: ε₀ = 0.95 + 0.01×LAI  (LAI ≤ 3) → RL↑ uchun
      ε₀ = 0.98             (LAI > 3)

F.22: εNB = 0.97 + 0.0033×LAI (LAI ≤ 3) → Ts uchun (L2 da kerak emas)
      εNB = 0.98              (LAI > 3)

Suv/qor (NDVI ≤ 0): ε₀ = εNB = 0.985
"""
import ee


class EmissivityCalculator:
    """Yer yuzasi emissiviteti."""
    
    def calc_emissivity(self, lai, ndvi):
        """
        ε₀ va εNB hisoblash.
        
        Returns:
            tuple: (e0, enb) — ee.Image juftligi
        """
        e0 = self._calc_e0(lai, ndvi)
        enb = self._calc_enb(lai, ndvi)
        return e0, enb
    
    def _calc_e0(self, lai, ndvi):
        """F.17: Keng-band emissivitet — RL↑ hisoblash uchun."""
        # LAI ≤ 3: ε₀ = 0.95 + 0.01 × LAI
        e0 = lai.multiply(0.01).add(0.95)
        
        # LAI > 3: ε₀ = 0.98
        e0 = e0.where(lai.gt(3), 0.98)
        
        # Suv/qor (NDVI ≤ 0): ε₀ = 0.985
        e0 = e0.where(ndvi.lte(0), 0.985)
        
        return e0.rename('e0')
    
    def _calc_enb(self, lai, ndvi):
        """F.22: Tor-band emissivitet — L2 Ts da kerak emas, lekin boshqa joyda foydali."""
        # LAI ≤ 3: εNB = 0.97 + 0.0033 × LAI
        enb = lai.multiply(0.0033).add(0.97)
        
        # LAI > 3: εNB = 0.98
        enb = enb.where(lai.gt(3), 0.98)
        
        # Suv/qor
        enb = enb.where(ndvi.lte(0), 0.985)
        
        return enb.rename('enb')
