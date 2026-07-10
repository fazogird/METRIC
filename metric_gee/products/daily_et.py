"""
Module 11: Kunlik ET.

F.54: ETrF = ETinst / ETr
F.55: ET24  = Crad × ETrF × ETr24
F.56: Crad  = (Rso_inst_H / Rso_inst_P) × (Rso_24_P / Rso_24_H)
F.57: Rso_24 = ∫₀²⁴ Rso_i dt  (24-soatlik aniq osmon radiatsiyasi)

DIQQAT — tau_sw va d² Crad nisbatida bekor bo'ladi:
  Rso = Gsc × cos(θ) × τsw / d²
  → Crad = cos(θhor)/cos(θrel) × Σcos(θrel_i)/Σcos(θhor_i)
  Faqat cos(θ) nisbatlari kerak.

F.7: cos(θrel) — qiyalik piksel uchun:
  cos(θrel) = sinδ·sinφ·coss − sinδ·cosφ·sins·cosγ
            + cosδ·cosφ·coss·cosω + cosδ·sinφ·sins·cosγ·cosω
            + cosδ·sins·sinγ·sinω

F.8: cos(θhor) — tekis yer uchun:
  cos(θhor) = sinδ·sinφ + cosδ·cosφ·cosω

Aspect konversiyasi (GEE → METRIC γ):
  GEE:    0°=N, 90°=E, 180°=S, 270°=W  (soat yo'nalishida)
  METRIC: γ=0→Janub, γ=±π→Shimol, γ=+π/2→G'arb, γ=-π/2→Sharq
  cos(γ) = -cos(aspect_gee_rad)
  sin(γ) =  sin(aspect_gee_rad)
"""
import math
import ee


class DailyET:
    """Kunlik bug'lanish hisoblash — F.54-57."""

    # ================================================================
    # F.54 — ETrF
    # ================================================================

    def calc_etrf(self, et_inst, etr):
        """F.54: ETrF = ETinst / ETr — crop coefficient analog."""
        etr_matched = (
            etr
            .resample('bilinear')
            .reproject(crs=et_inst.projection())
            .max(0.01)
        )
        return (et_inst
                .divide(etr_matched.max(0.01))
                .clamp(0, 1.05)
                .rename('ETrF'))

    # ================================================================
    # F.55 — ET24
    # ================================================================

    def calc_et24(self, etrf, etr_24, crad=None):
        """
        F.55: ET24 = Crad × ETrF × ETr24  (mm/day)

        Args:
            etrf   : ee.Image — ETrF (dimensionless)
            etr_24 : ee.Image — ETr 24-soatlik (mm/day)
            crad   : ee.Image — qiyalik tuzatma (F.56); None → 1.0 (tekis yer)
        """
        if crad is None:
            crad = ee.Image.constant(1.0)
        return (crad
                .multiply(etrf)
                .multiply(etr_24)
                .rename('ET24'))

    # ================================================================
    # F.56 — Crad
    # ================================================================

    def calc_crad(self, dem, solar):
        """
        F.56: Crad = (Rso_inst_H / Rso_inst_P) × (Rso_24_P / Rso_24_H)

        Birinchi nisbat — lahzali (image overpass vaqti):
          tau_sw va d² bekor bo'ladi → cos(θhor) / cos(θrel)

        Ikkinchi nisbat — kunlik integral (F.57):
          Σcos(θrel_i) / Σcos(θhor_i)  i ∈ [0..23] soat

        Args:
            dem   : ee.Image — 'slope_rad', 'aspect_rad' bandlari (DEM loader dan)
            solar : dict     — 'delta'(ee.Number), 'omega'(ee.Number),
                               'cos_theta_hor'(ee.Number)

        Returns:
            ee.Image — 'Crad', clamp [0.5, 1.5]
        """
        # ── Lahzali nisbat ─────────────────────────────────────────
        # cos(θhor) — ee.Number → ee.Image
        cos_h_inst = solar['cos_theta_hor'].max(0.001)

        # F.7 — cos(θrel) image overpass vaqtida
        cos_p_inst = self._calc_cos_theta_slope(dem, solar).max(0.001)

        rso_inst_ratio = cos_h_inst.divide(cos_p_inst)    # Rso_inst_H / Rso_inst_P

        # ── F.57: 24-soatlik yig'indi nisbati ──────────────────────
        rso_24_h, rso_24_p = self._calc_rso_24(dem, solar)
        rso_24_ratio = rso_24_p.divide(rso_24_h.max(0.001))  # Rso_24_P / Rso_24_H

        # ── F.56: Crad ─────────────────────────────────────────────
        crad = rso_inst_ratio.multiply(rso_24_ratio)
        return crad.clamp(0.5, 1.5).rename('Crad')

    # ================================================================
    # F.7 — cos(θrel): qiyalik uchun burchak
    # ================================================================

    def _calc_cos_theta_slope(self, dem, solar):
        """
        F.7: cos(θrel) — qiyalik piksel uchun quyosh tushish burchagi.

        5 ta had (Allen 2007):
          T1 =  sinδ · sinφ · cos(s)
          T2 = -sinδ · cosφ · sin(s) · cos(γ)
          T3 =  cosδ · cosφ · cos(s) · cos(ω)
          T4 =  cosδ · sinφ · sin(s) · cos(γ) · cos(ω)
          T5 =  cosδ · sin(s) · sin(γ) · sin(ω)

        Args:
            dem   : ee.Image — 'slope_rad', 'aspect_rad' (DEM loader dan)
            solar : dict     — 'delta'(rad, ee.Number), 'omega'(rad, ee.Number)

        Returns:
            ee.Image — cos(θrel), manfiy bo'lishi mumkin (soya)
        """
        delta = solar['delta']   # ee.Number — quyosh deklinatsiyasi (rad)
        omega = solar['omega']   # ee.Number — soat burchagi (rad)

        lat_rad = ee.Image.pixelLonLat().select('latitude').multiply(math.pi / 180)

        slope_s = dem.select('slope_rad')       # s (rad)
        asp_rad = dem.select('aspect_rad')       # GEE convention (rad)

        # GEE → METRIC γ konversiyasi:
        #   cos(γ) = -cos(aspect_gee)
        #   sin(γ) =  sin(aspect_gee)
        cos_gamma = asp_rad.cos().multiply(-1)
        sin_gamma = asp_rad.sin()

        # Hisob
        sin_d = delta.sin()
        cos_d = delta.cos()
        sin_p = lat_rad.sin()
        cos_p = lat_rad.cos()
        cos_s = slope_s.cos()
        sin_s = slope_s.sin()
        cos_w = omega.cos()
        sin_w = omega.sin()

        # ⚠️ GEE Python: Number.multiply(Image) → Number!
        #               Image.multiply(Number) → Image ← to'g'ri
        t1 = sin_p.multiply(sin_d).multiply(cos_s)
        t2 = cos_p.multiply(sin_d).multiply(sin_s).multiply(cos_gamma).multiply(-1)
        t3 = cos_p.multiply(cos_d).multiply(cos_s).multiply(cos_w)
        t4 = sin_p.multiply(cos_d).multiply(sin_s).multiply(cos_gamma).multiply(cos_w)
        t5 = sin_s.multiply(cos_d).multiply(sin_gamma).multiply(sin_w)

        return t1.add(t2).add(t3).add(t4).add(t5).rename('cos_theta_rel')

    # ================================================================
    # F.57 — Rso_24: 24-soatlik aniq osmon radiatsiyasi (nisbatlar)
    # ================================================================

    def _calc_rso_24(self, dem, solar):
        """
        F.57: Rso_24 = Σᵢ Rso_i  (i = 0..23 soat)

        Tau_sw va d² — o'zgarmas, Crad nisbatida bekor bo'ladi.
        Shuning uchun faqat cos(θ) yig'indisi hisoblanadi.

        Har soat uchun omega_i = (i + 0.5 − 12) × π/12  rad
        (soat markazi, −11.5 soat ... +11.5 soat)

        Returns:
            rso_24_h : ee.Image — tekis yer Σcos(θhor_i), faqat kun (>0)
            rso_24_p : ee.Image — qiyalik  Σcos(θrel_i), faqat kun (>0)
        """
        delta = solar['delta']   # ee.Number

        lat_rad = ee.Image.pixelLonLat().select('latitude').multiply(math.pi / 180)

        slope_s = dem.select('slope_rad')
        asp_rad = dem.select('aspect_rad')
        cos_gamma = asp_rad.cos().multiply(-1)
        sin_gamma = asp_rad.sin()

        sin_d = delta.sin()
        cos_d = delta.cos()
        sin_p = lat_rad.sin()
        cos_p = lat_rad.cos()
        cos_s = slope_s.cos()
        sin_s = slope_s.sin()

        rso_24_h = ee.Image.constant(0)
        rso_24_p = ee.Image.constant(0)

        for h in range(24):
            # Soat markazidagi soat burchagi
            omega_h = ee.Number((h + 0.5 - 12) * math.pi / 12)
            cos_w = omega_h.cos()
            sin_w = omega_h.sin()

            # F.8: tekis yer cos(θhor)
            # ⚠️ Image CHAPDA: Image.multiply(Number) → Image
            cos_h = (sin_p.multiply(sin_d)
                     .add(cos_p.multiply(cos_d).multiply(cos_w)))

            # F.7: qiyalik cos(θrel)
            t1 = sin_p.multiply(sin_d).multiply(cos_s)
            t2 = cos_p.multiply(sin_d).multiply(sin_s).multiply(cos_gamma).multiply(-1)
            t3 = cos_p.multiply(cos_d).multiply(cos_s).multiply(cos_w)
            t4 = sin_p.multiply(cos_d).multiply(sin_s).multiply(cos_gamma).multiply(cos_w)
            t5 = sin_s.multiply(cos_d).multiply(sin_gamma).multiply(sin_w)
            cos_p_h = t1.add(t2).add(t3).add(t4).add(t5)

            # Faqat kun vaqti: cos > 0 (quyosh ufqdan yuqori)
            rso_24_h = rso_24_h.add(cos_h.max(0))
            rso_24_p = rso_24_p.add(cos_p_h.max(0))

        return rso_24_h, rso_24_p
