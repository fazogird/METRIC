"""
CIMEC Anchor Selector — LULC asosida cold/hot piksel tanlash.

CASCADE FALLBACK (tartib bilan):
  L1 CIMEC   : LULC + NDVI P95 + Ts P20   /  NDVI P10 + Ts P80
  L2 Plan A  : LULC + NDVI P80 + Ts P30   /  NDVI P20 + Ts P70
  L3 Plan B  : LULC + NDVI P50 + Ts P40   /  NDVI P30 + Ts P60
  L4 pySEBAL : Butun area, LULC yo'q (oxirgi fallback)

ESA WorldCover v100 (2020) kodlari:
  10  Trees             → yo'q
  20  Shrubland         → yo'q
  30  Grassland         → yo'q
  40  Cropland          → COLD (sug'oriladigan ekin)
  50  Built-up          → HOT  (shahar/yo'l)
  60  Bare/sparse veg   → HOT  (yalang'och tuproq / cho'l)
  70  Snow and ice      → yo'q
  80  Permanent water   → yo'q
  90  Herbaceous wetland→ yo'q
  95  Mangroves         → yo'q
  100 Moss/lichen       → yo'q

Markaziy Osiyo uchun ham ishlaydi:
  cold: [40]        → sug'oriladigan dalalar
  hot:  [60, 50]    → cho'l, yalang'och tuproq, shahar
"""
import ee
from typing import List, Tuple


# ─── ESA WorldCover ───────────────────────────────────────────────
ESA_WORLDCOVER = "ESA/WorldCover/v100/2020"

# Default LULC kodlar
DEFAULT_COLD_LC = [40]          # Cropland
DEFAULT_HOT_LC  = [60, 50]      # Bare + Built-up


class AnchorSelector:
    """
    CIMEC anchor piksel tanlash — LULC + cascade fallback.

    Ishlatish:
        selector = AnchorSelector(
            cold_lc=[40],
            hot_lc=[60, 50],
            lapse_rate=6.5
        )
        result = selector.select(ndvi, ts_dem, geometry)
        cold_mask = result['cold_mask']
        hot_mask  = result['hot_mask']
        plan_used = result['plan']
    """

    def __init__(
        self,
        cold_lc: List[int] = None,
        hot_lc:  List[int] = None,
        lapse_rate: float = 6.5,
        verbose: bool = True,
    ):
        self.cold_lc    = cold_lc or DEFAULT_COLD_LC
        self.hot_lc     = hot_lc  or DEFAULT_HOT_LC
        self.lapse_rate = lapse_rate   # K/km (DEM korreksiya uchun)
        self.verbose    = verbose

    # ================================================================
    # ASOSIY METOD
    # ================================================================

    def select(self, ndvi, ts_dem, geometry) -> dict:
        """
        Cold/hot piksel maskalarini qaytaradi.

        Args:
            ndvi     : ee.Image — NDVI (band nomi 'NDVI')
            ts_dem   : ee.Image — lapse corrected Ts (band nomi 'TsDEM')
            geometry : ee.Geometry — qidiruv hududi

        Returns:
            dict:
                'cold_mask' : ee.Image
                'hot_mask'  : ee.Image
                'cold_ts'   : ee.Number — cold anchor TsDEM o'rtacha
                'hot_ts'    : ee.Number — hot anchor TsDEM o'rtacha
                'plan'      : str       — ishlatilgan metod nomi
                'cold_ndvi' : ee.Number
                'hot_ndvi'  : ee.Number
        """
        # LULC xaritasi (ESA WorldCover 2020)
        lulc = ee.Image(ESA_WORLDCOVER).select('Map').clip(geometry)

        cold_lulc_mask = self._build_lc_mask(lulc, self.cold_lc)
        hot_lulc_mask  = self._build_lc_mask(lulc, self.hot_lc)

        lulc_pixel_count = lulc.reduceRegion(
            ee.Reducer.count(), geometry, 100, maxPixels=1e9
        ).values().get(0)
        lulc_count = ee.Number(lulc_pixel_count).getInfo()
        self._log(f"  LULC piksel soni: {lulc_count}")
        self._log(f"  Cold LC kodlar: {self.cold_lc}   Hot LC kodlar: {self.hot_lc}")

        # ── CASCADE ──────────────────────────────────────────────
        plans = [
            ("L1-CIMEC",   self._cimec_cold,    self._cimec_hot,    cold_lulc_mask, hot_lulc_mask),
            ("L2-PlanA",   self._plan_a_cold,   self._plan_a_hot,   cold_lulc_mask, hot_lulc_mask),
            ("L3-PlanB",   self._plan_b_cold,   self._plan_b_hot,   cold_lulc_mask, hot_lulc_mask),
            ("L4-pySEBAL", self._pysebal_cold,  self._pysebal_hot,  None,           None),
        ]

        for name, cold_fn, hot_fn, c_lc, h_lc in plans:
            try:
                cold_mask = cold_fn(ndvi, ts_dem, geometry, c_lc)
                hot_mask  = hot_fn( ndvi, ts_dem, geometry, h_lc)

                # Piksel sonini tekshirish
                n_cold = self._pixel_count(cold_mask, geometry)
                n_hot  = self._pixel_count(hot_mask,  geometry)

                self._log(f"  {name}: cold={n_cold}  hot={n_hot}")

                if n_cold > 0 and n_hot > 0:
                    self._log(f"  ✅ Anchor tanlandi: {name}")
                    cold_ts, cold_ndvi = self._anchor_stats(ts_dem, ndvi, cold_mask, geometry)
                    hot_ts,  hot_ndvi  = self._anchor_stats(ts_dem, ndvi, hot_mask,  geometry)
                    return {
                        'cold_mask': cold_mask,
                        'hot_mask':  hot_mask,
                        'cold_ts':   cold_ts,
                        'hot_ts':    hot_ts,
                        'cold_ndvi': cold_ndvi,
                        'hot_ndvi':  hot_ndvi,
                        'plan':      name,
                    }
                else:
                    self._log(f"  ⏩ {name}: piksel yo'q → keyingi")

            except Exception as e:
                self._log(f"  ⚠️ {name} xato: {str(e)[:80]} → keyingi")

        # Hech qaysi ishlmasa — xato
        raise RuntimeError("Hech qaysi anchor metod cold/hot piksel topa olmadi!")

    # ================================================================
    # LEVEL 1 — CIMEC (Allen 2007 standart)
    # ================================================================

    def _cimec_cold(self, ndvi, ts_dem, geometry, lc_mask):
        """Top 5% NDVI → coldest 20% TsDEM (LULC ichida)."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p95 = self._percentile(ndvi_m, 'NDVI', 95, geometry)
        top = ndvi_m.gte(p95)
        p20 = self._percentile(ts_m.updateMask(top), 'TsDEM', 20, geometry)
        return top.And(ts_m.lte(p20)).selfMask()

    def _cimec_hot(self, ndvi, ts_dem, geometry, lc_mask):
        """Lowest 10% NDVI → hottest 20% TsDEM (LULC ichida)."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p10 = self._percentile(ndvi_m, 'NDVI', 10, geometry)
        low = ndvi_m.lte(p10)
        p80 = self._percentile(ts_m.updateMask(low), 'TsDEM', 80, geometry)
        return low.And(ts_m.gte(p80)).selfMask()

    # ================================================================
    # LEVEL 2 — Plan A (kengroq threshold)
    # ================================================================

    def _plan_a_cold(self, ndvi, ts_dem, geometry, lc_mask):
        """Top 20% NDVI → coldest 30% Ts."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p80 = self._percentile(ndvi_m, 'NDVI', 80, geometry)
        top = ndvi_m.gte(p80)
        p30 = self._percentile(ts_m.updateMask(top), 'TsDEM', 30, geometry)
        return top.And(ts_m.lte(p30)).selfMask()

    def _plan_a_hot(self, ndvi, ts_dem, geometry, lc_mask):
        """Lowest 20% NDVI → hottest 30% Ts."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p20 = self._percentile(ndvi_m, 'NDVI', 20, geometry)
        low = ndvi_m.lte(p20)
        p70 = self._percentile(ts_m.updateMask(low), 'TsDEM', 70, geometry)
        return low.And(ts_m.gte(p70)).selfMask()

    # ================================================================
    # LEVEL 3 — Plan B (yanada kengroq)
    # ================================================================

    def _plan_b_cold(self, ndvi, ts_dem, geometry, lc_mask):
        """Top 50% NDVI → coldest 40% Ts."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p50 = self._percentile(ndvi_m, 'NDVI', 50, geometry)
        top = ndvi_m.gte(p50)
        p40 = self._percentile(ts_m.updateMask(top), 'TsDEM', 40, geometry)
        return top.And(ts_m.lte(p40)).selfMask()

    def _plan_b_hot(self, ndvi, ts_dem, geometry, lc_mask):
        """Lowest 30% NDVI → hottest 40% Ts."""
        ndvi_m = ndvi.updateMask(lc_mask) if lc_mask else ndvi
        ts_m   = ts_dem.updateMask(lc_mask) if lc_mask else ts_dem

        p30 = self._percentile(ndvi_m, 'NDVI', 30, geometry)
        low = ndvi_m.lte(p30)
        p60 = self._percentile(ts_m.updateMask(low), 'TsDEM', 60, geometry)
        return low.And(ts_m.gte(p60)).selfMask()

    # ================================================================
    # LEVEL 4 — pySEBAL fallback (LULC yo'q, butun area)
    # ================================================================

    def _pysebal_cold(self, ndvi, ts_dem, geometry, lc_mask=None):
        """Butun area: Top 5% NDVI → coldest 20% Ts."""
        p95 = self._percentile(ndvi, 'NDVI', 95, geometry)
        top = ndvi.gte(p95)
        p20 = self._percentile(ts_dem.updateMask(top), 'TsDEM', 20, geometry)
        return top.And(ts_dem.lte(p20)).selfMask()

    def _pysebal_hot(self, ndvi, ts_dem, geometry, lc_mask=None):
        """Butun area: Lowest 10% NDVI → hottest 20% Ts."""
        p10 = self._percentile(ndvi, 'NDVI', 10, geometry)
        low = ndvi.lte(p10)
        p80 = self._percentile(ts_dem.updateMask(low), 'TsDEM', 80, geometry)
        return low.And(ts_dem.gte(p80)).selfMask()

    # ================================================================
    # YORDAMCHI METODLAR
    # ================================================================

    @staticmethod
    def _build_lc_mask(lulc_img, codes: List[int]) -> ee.Image:
        """LULC kodlardan binary maska yasash."""
        mask = lulc_img.eq(codes[0])
        for c in codes[1:]:
            mask = mask.Or(lulc_img.eq(c))
        return mask.selfMask()

    @staticmethod
    def _percentile(img, band: str, pct: int, geometry) -> ee.Number:
        """
        GEE da percentile olish.
        Bir yoki ko'p percentile — kalit nomini auto topadi.
        """
        res = img.reduceRegion(
            ee.Reducer.percentile([pct]), geometry, 100, maxPixels=1e9
        ).getInfo()
        # Kalit: 'NDVI_p95' yoki 'NDVI' (bitta percentile)
        for key in [f'{band}_p{pct}', band]:
            if key in res and res[key] is not None:
                return ee.Number(res[key])
        raise ValueError(f"_percentile: '{band}' p{pct} topilmadi. Keys={list(res.keys())}")

    @staticmethod
    def _pixel_count(mask: ee.Image, geometry) -> int:
        """Maskadagi piksel sonini qaytaradi."""
        res = mask.reduceRegion(
            ee.Reducer.count(), geometry, 100, maxPixels=1e9
        ).values().get(0)
        return int(ee.Number(res).getInfo() or 0)

    @staticmethod
    def _anchor_stats(ts_dem, ndvi, mask, geometry) -> Tuple[ee.Number, ee.Number]:
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
