"""
CIMEC Anchor Selector — Allen et al. (2013) to'liq implementatsiya.

JARAYON (maqola tartibida):

HO'L (Cold) piksel:
  1. Top 5% NDVI (LULC filtrlangan)               [L1 CIMEC]
  2. Eng sovuq 20% TsDEM
  3. TsDEM ± 0.2 K toleransiya filtri             [Allen 2013]
  4. Albedo filtri: |α - α_limit| ≤ 0.02  [F.7]  [Allen 2013]

QURUQ (Hot) piksel:
  1. Eng past 10% NDVI (LULC filtrlangan)          [L1 CIMEC]
  2. Eng issiq 20% TsDEM
  3. TsDEM ± 0.2 K toleransiya filtri             [Allen 2013]
  4. Tfac tuzatmasi [F.8]:
       Tfac = max(0, 2.6 - 13 × P₆₀/ETr₆₀)
       hot_ts_corrected = hot_ts - Tfac

F.7:  α_limit = 0.001343·β + 0.3281·exp(-0.0188·β)
F.8:  Tfac    = max(0, 2.6 - 13·P₆₀/ETr₆₀)  — faqat P/ETr < 0.2 bo'lganda

ETrF tayinlash (pipeline.py da):
  F.4a: ETrFcold = 1.54·NDVI₁ - 0.10  (NDVI₁ < 0.75)
  F.4b: ETrFcold = 1.05               (NDVI₁ ≥ 0.75)
  F.5:  fc = (NDVIhot - NDVIbare) / (NDVIfull - NDVIbare)
  F.6:  ETrFhot = fc·ETrFcold + (1-fc)·ETrFbare

CASCADE FALLBACK:
  L1 CIMEC   : LULC + P95 NDVI + P20 Ts + ±0.2K + albedo [F.7] + Tfac [F.8]
  L2 Plan A  : LULC + P80 NDVI + P30 Ts  (strict filtrlar yo'q)
  L3 Plan B  : LULC + P50 NDVI + P40 Ts  (strict filtrlar yo'q)
  L4 pySEBAL : Butun area, LULC yo'q      (strict filtrlar yo'q)

ESA WorldCover v100 (2020) kodlari:
  40  Cropland          → COLD (sug'oriladigan ekin)
  50  Built-up          → HOT  (shahar/yo'l)
  60  Bare/sparse veg   → HOT  (yalang'och tuproq / cho'l)
"""
import ee
from typing import List, Optional, Tuple


# ── ESA WorldCover ──────────────────────────────────────────────────
ESA_WORLDCOVER = "ESA/WorldCover/v100/2020"
DEFAULT_COLD_LC = [40]       # Cropland
DEFAULT_HOT_LC  = [60, 50]   # Bare + Built-up

# Allen 2013 toleransiyalar (constants.py da ham mavjud)
TS_TOLERANCE  = 0.2    # K  — TsDEM toleransiya filtri
ALBEDO_TOL    = 0.02   # —  — Albedo toleransiya filtri


class AnchorSelector:
    """
    CIMEC anchor piksel tanlash — Allen et al. (2013) formulalar bilan.

    Ishlatish:
        selector = AnchorSelector(cold_lc=[40], hot_lc=[60,50])
        result = selector.select(ndvi, ts_dem, geometry,
                                 albedo=albedo_img,
                                 sun_elev=sun_elev_val,
                                 precip_60_val=p60,
                                 etr_60_val=e60)
        cold_mask        = result['cold_mask']
        hot_ts_corrected = result['hot_ts_corrected']   # F.8 tuzatilgan
        cold_ndvi        = result['cold_ndvi']           # F.4a uchun
    """

    def __init__(
        self,
        cold_lc:    List[int] = None,
        hot_lc:     List[int] = None,
        lapse_rate: float = 6.5,
        verbose:    bool  = True,
    ):
        self.cold_lc    = cold_lc or DEFAULT_COLD_LC
        self.hot_lc     = hot_lc  or DEFAULT_HOT_LC
        self.lapse_rate = lapse_rate
        self.verbose    = verbose

    # ================================================================
    # ASOSIY METOD
    # ================================================================

    def select(
        self,
        ndvi,
        ts_dem,
        geometry,
        albedo=None,           # F.7: sovuq piksel albedo filtri uchun (ee.Image)
        sun_elev=None,         # F.7: β quyosh balandligi (daraja, float yoki ee.Number)
        precip_60_val=None,    # F.8: 60-kunlik yog'in (mm, float)
        etr_60_val=None,       # F.8: 60-kunlik ETr   (mm, float)
    ) -> dict:
        """
        Maqola tartibida cold/hot anchor piksel maskalarini tanlaydi.

        Returns:
            dict:
              cold_mask        : ee.Image
              hot_mask         : ee.Image
              cold_ts          : ee.Number — sovuq anchor TsDEM o'rtacha [K]
              hot_ts           : ee.Number — issiq anchor fizik TsDEM [K]  (λ uchun)
              hot_ts_corrected : ee.Number — hot_ts - Tfac [K]  (regression uchun, F.8)
              tfac             : ee.Number — Tfac qiymati [K]
              cold_ndvi        : ee.Number — F.4a uchun
              hot_ndvi         : ee.Number — F.5/F.6 uchun
              plan             : str       — ishlatilgan cascade nomi
        """
        # LULC xaritasi
        lulc       = ee.Image(ESA_WORLDCOVER).select('Map').clip(geometry)
        cold_lulc  = self._build_lc_mask(lulc, self.cold_lc)
        hot_lulc   = self._build_lc_mask(lulc, self.hot_lc)

        self._log(f"  Cold LC: {self.cold_lc}   Hot LC: {self.hot_lc}")

        # ── F.7: Albedo chegarasi ─────────────────────────────────
        alpha_limit = None
        if albedo is not None and sun_elev is not None:
            beta        = float(ee.Number(sun_elev).getInfo()
                                if isinstance(sun_elev, ee.ComputedObject) else sun_elev)
            alpha_limit = self._calc_alpha_limit(beta)
            self._log(f"  F.7 α_limit={alpha_limit:.4f}  (β={beta:.1f}°)")

        # ── F.8: Tfac ─────────────────────────────────────────────
        tfac = ee.Number(0)
        if precip_60_val is not None and etr_60_val is not None:
            p60 = float(precip_60_val)
            e60 = float(etr_60_val)
            tfac_val = self._calc_tfac_py(p60, e60)
            tfac     = ee.Number(tfac_val)
            self._log(f"  F.8 Tfac={tfac_val:.2f} K  "
                      f"(P₆₀={p60:.1f} mm  ETr₆₀={e60:.1f} mm  P/ETr={p60/max(e60,1):.3f})")

        # ── CASCADE ───────────────────────────────────────────────
        # (name, cold_lc, hot_lc, ndvi_c_pct, ts_c_pct, ndvi_h_pct, ts_h_pct, strict)
        levels = [
            ("L1-CIMEC",   cold_lulc, hot_lulc, 95, 20, 10, 80, True),
            ("L2-PlanA",   cold_lulc, hot_lulc, 80, 30, 20, 70, False),
            ("L3-PlanB",   cold_lulc, hot_lulc, 50, 40, 30, 60, False),
            ("L4-pySEBAL", None,      None,     95, 20, 10, 80, False),
        ]

        for name, c_lc, h_lc, ndvi_c, ts_c, ndvi_h, ts_h, strict in levels:
            try:
                # Cold piksel tanlash
                cold_mask = self._select_cold(
                    ndvi, ts_dem, geometry, c_lc,
                    ndvi_pct  = ndvi_c,
                    ts_pct    = ts_c,
                    albedo    = albedo      if strict else None,
                    alpha_lim = alpha_limit if strict else None,
                    strict    = strict,
                )
                # Hot piksel tanlash
                hot_mask, hot_ts_phys, hot_ts_corr = self._select_hot(
                    ndvi, ts_dem, geometry, h_lc,
                    ndvi_pct = ndvi_h,
                    ts_pct   = ts_h,
                    tfac     = tfac if strict else ee.Number(0),
                    strict   = strict,
                )

                n_cold = self._pixel_count(cold_mask, geometry)
                n_hot  = self._pixel_count(hot_mask,  geometry)
                self._log(f"  {name}: cold={n_cold}  hot={n_hot}")

                if n_cold > 0 and n_hot > 0:
                    self._log(f"  ✅ Anchor tanlandi: {name}")

                    cold_ts, cold_ndvi = self._anchor_stats(
                        ts_dem, ndvi, cold_mask, geometry)
                    _,       hot_ndvi  = self._anchor_stats(
                        ts_dem, ndvi, hot_mask,  geometry)

                    return {
                        'cold_mask':        cold_mask,
                        'hot_mask':         hot_mask,
                        'cold_ts':          cold_ts,
                        'hot_ts':           hot_ts_phys,    # fizik (λ uchun)
                        'hot_ts_corrected': hot_ts_corr,    # F.8 tuzatilgan (regression uchun)
                        'tfac':             tfac if strict else ee.Number(0),
                        'cold_ndvi':        cold_ndvi,
                        'hot_ndvi':         hot_ndvi,
                        'plan':             name,
                    }
                else:
                    self._log(f"  ⏩ {name}: piksel yo'q → keyingi")

            except Exception as e:
                self._log(f"  ⚠️ {name}: {str(e)[:80]} → keyingi")

        raise RuntimeError("Hech qaysi anchor metod cold/hot piksel topa olmadi!")

    # ================================================================
    # COLD PIKSEL TANLASH
    # ================================================================

    def _select_cold(
        self,
        ndvi, ts_dem, geometry, lc_mask,
        ndvi_pct: int = 95,
        ts_pct:   int = 20,
        albedo    = None,
        alpha_lim = None,
        strict:   bool = False,
    ):
        """
        Maqola tartibida cold piksel tanlash:
          1. Top X% NDVI (LULC filtrlangan)
          2. Eng sovuq Y% TsDEM
          3. TsDEM ± 0.2 K toleransiya filtri  [strict=True]
          4. Albedo filtri |α - α_limit| ≤ 0.02  [F.7, strict=True]

        Returns:
            ee.Image — cold_mask (selfMasked)
        """
        ndvi_m = ndvi.updateMask(lc_mask)  if lc_mask is not None else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask is not None else ts_dem

        # 1. Top X% NDVI
        p_ndvi = self._percentile(ndvi_m, 'NDVI', ndvi_pct, geometry)
        top    = ndvi_m.gte(p_ndvi)

        # 2. Eng sovuq Y% TsDEM
        p_ts      = self._percentile(ts_m.updateMask(top), 'TsDEM', ts_pct, geometry)
        cold_mask = top.And(ts_m.lte(p_ts))

        if strict:
            # 3. TsDEM ± 0.2 K toleransiya filtri
            ts_mean_val = (ts_m.updateMask(cold_mask)
                          .reduceRegion(ee.Reducer.mean(), geometry, 100, maxPixels=1e9)
                          .getInfo().get('TsDEM'))
            if ts_mean_val is not None:
                cold_mask = cold_mask.And(
                    ts_m.subtract(ts_mean_val).abs().lte(TS_TOLERANCE)
                )

            # 4. F.7 Albedo filtri
            if albedo is not None and alpha_lim is not None:
                albedo_filter = albedo.subtract(alpha_lim).abs().lte(ALBEDO_TOL)
                cold_mask     = cold_mask.And(albedo_filter)

        return cold_mask.selfMask()

    # ================================================================
    # HOT PIKSEL TANLASH
    # ================================================================

    def _select_hot(
        self,
        ndvi, ts_dem, geometry, lc_mask,
        ndvi_pct: int = 10,
        ts_pct:   int = 80,
        tfac      = None,
        strict:   bool = False,
    ) -> Tuple:
        """
        Maqola tartibida hot piksel tanlash:
          1. Eng past X% NDVI (LULC filtrlangan)
          2. Eng issiq Y% TsDEM
          3. TsDEM ± 0.2 K toleransiya filtri  [strict=True]
          4. Tfac tuzatmasi [F.8, strict=True]:
               hot_ts_corrected = hot_ts_mean - Tfac

        Returns:
            (hot_mask, hot_ts_physical, hot_ts_corrected)
              hot_ts_physical  — o'lchangan fizik TsDEM (λ hisoblash uchun)
              hot_ts_corrected — Tfac bilan tuzatilgan TsDEM (regression uchun)
        """
        ndvi_m = ndvi.updateMask(lc_mask)  if lc_mask is not None else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask is not None else ts_dem

        # 1. Eng past X% NDVI
        p_ndvi   = self._percentile(ndvi_m, 'NDVI', ndvi_pct, geometry)
        low      = ndvi_m.lte(p_ndvi)

        # 2. Eng issiq Y% TsDEM
        p_ts    = self._percentile(ts_m.updateMask(low), 'TsDEM', ts_pct, geometry)
        hot_mask = low.And(ts_m.gte(p_ts))

        # Boshlang'ich o'rtacha fizik TsDEM
        ts_mean_val = (ts_m.updateMask(hot_mask)
                      .reduceRegion(ee.Reducer.mean(), geometry, 100, maxPixels=1e9)
                      .getInfo().get('TsDEM')) or 310.0
        hot_ts_physical = ee.Number(ts_mean_val)

        if strict:
            # 3. TsDEM ± 0.2 K toleransiya filtri
            hot_mask = hot_mask.And(
                ts_m.subtract(hot_ts_physical).abs().lte(TS_TOLERANCE)
            )
            # Yangilangan o'rtacha (faqat toleransiya doirasida)
            ts_mean_val2 = (ts_m.updateMask(hot_mask)
                           .reduceRegion(ee.Reducer.mean(), geometry, 100, maxPixels=1e9)
                           .getInfo().get('TsDEM'))
            if ts_mean_val2 is not None:
                hot_ts_physical = ee.Number(ts_mean_val2)

        # 4. F.8: Tfac tuzatmasi (regression uchun)
        _tfac           = tfac if tfac is not None else ee.Number(0)
        hot_ts_corrected = hot_ts_physical.subtract(_tfac)

        return hot_mask.selfMask(), hot_ts_physical, hot_ts_corrected

    # ================================================================
    # F.7 — Albedo chegarasi  (Python skalyar, cold cluster uchun)
    # ================================================================

    @staticmethod
    def _calc_alpha_limit(beta_deg: float) -> float:
        """
        F.7: α_limit = 0.001343·β + 0.3281·exp(-0.0188·β)

        Args:
            beta_deg: float — β, quyosh balandligi (daraja)
        Returns:
            float — α_limit
        """
        import math
        return 0.001343 * beta_deg + 0.3281 * math.exp(-0.0188 * beta_deg)

    # ================================================================
    # F.8 — Tfac  (Python skalyar)
    # ================================================================

    @staticmethod
    def _calc_tfac_py(precip_60: float, etr_60: float) -> float:
        """
        F.8: Tfac = max(0, 2.6 - 13 × P₆₀/ETr₆₀)

        Faqat P/ETr < 0.2 bo'lsa tuzatma kerak (qurg'oqchilik sharoit).

        Args:
            precip_60: float — 60 kunlik yog'in (mm)
            etr_60   : float — 60 kunlik ETr (mm)
        Returns:
            float — Tfac [K], ≥ 0
        """
        etr_safe = max(etr_60, 1.0)
        p_ratio  = precip_60 / etr_safe
        if p_ratio >= 0.2:
            return 0.0           # Ho'l sharoit — tuzatma kerak emas
        tfac = max(0.0, 2.6 - 13.0 * p_ratio)
        return tfac

    # ================================================================
    # YORDAMCHI METODLAR
    # ================================================================

    @staticmethod
    def _build_lc_mask(lulc_img, codes: List[int]) -> 'ee.Image':
        """LULC kodlardan binary maska yasash."""
        mask = lulc_img.eq(codes[0])
        for c in codes[1:]:
            mask = mask.Or(lulc_img.eq(c))
        return mask.selfMask()

    @staticmethod
    def _percentile(img, band: str, pct: int, geometry) -> 'ee.Number':
        """GEE da percentile olish — kalit nomini auto topadi."""
        res = img.reduceRegion(
            ee.Reducer.percentile([pct]), geometry, 100, maxPixels=1e9
        ).getInfo()
        for key in [f'{band}_p{pct}', band]:
            if key in res and res[key] is not None:
                return ee.Number(res[key])
        raise ValueError(
            f"_percentile: '{band}' p{pct} topilmadi. Keys={list(res.keys())}")

    @staticmethod
    def _pixel_count(mask: 'ee.Image', geometry) -> int:
        """Maskadagi piksel sonini qaytaradi."""
        res = mask.reduceRegion(
            ee.Reducer.count(), geometry, 100, maxPixels=1e9
        ).values().get(0)
        return int(ee.Number(res).getInfo() or 0)

    @staticmethod
    def _anchor_stats(
        ts_dem, ndvi, mask, geometry
    ) -> Tuple['ee.Number', 'ee.Number']:
        """Anchor mask uchun TsDEM va NDVI o'rtachasi."""
        res = (ts_dem.addBands(ndvi)
               .updateMask(mask)
               .reduceRegion(ee.Reducer.mean(), geometry, 100, maxPixels=1e9)
               .getInfo())
        ts_val   = res.get('TsDEM') or res.get('TsDEM_mean') or 300.0
        ndvi_val = res.get('NDVI')  or res.get('NDVI_mean')  or 0.3
        return ee.Number(ts_val), ee.Number(ndvi_val)

    def _log(self, msg: str):
        if self.verbose:
            print(msg)
