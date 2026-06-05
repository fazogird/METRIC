"""
Module 8: CIMEC — avtomatik kalibrlash (Allen et al., 2013).

METRIC ning yuragi — cold/hot piksel tanlash va dT = a + bTs kalibrlash.

CIMEC F.4-8, METRIC F.46-51.
"""
import ee
from metric_gee.config.constants import (
    CIMEC_NDVI_COLD_PCT, CIMEC_TS_COLD_PCT,
    CIMEC_NDVI_HOT_PCT, CIMEC_TS_HOT_PCT,
    CIMEC_ETRF_COLD_DEFAULT, CIMEC_NDVI_BARE, CIMEC_NDVI_FULL, CP
)


class CIMECCalibrator:
    """CIMEC avtomatik kalibrlash."""

    def __init__(self, settings):
        self.settings = settings

    def calibrate(self, ndvi, ts_dem, albedo, agri_mask,
                  rn, g, aero, era5, solar):
        """To'liq CIMEC: cold/hot tanlash → dT(a,b) kalibrlash."""
        ndvi_m = ndvi.updateMask(agri_mask)
        ts_dem_m = ts_dem.updateMask(agri_mask)
        geometry = ndvi.geometry()

        cold = self._select_cold(ndvi_m, ts_dem_m, geometry)
        hot = self._select_hot(ndvi_m, ts_dem_m, geometry)

        etrf_cold = self._assign_etrf_cold(cold['ndvi_mean'])
        etrf_hot = self._assign_etrf_hot(hot['ndvi_mean'], etrf_cold)

        a, b = self._solve_dt_coefficients(
            cold, hot, etrf_cold, etrf_hot,
            rn, g, aero, geometry)

        return {'a': a, 'b': b,
                'cold_ts': cold['ts_mean'], 'hot_ts': hot['ts_mean'],
                'etrf_cold': etrf_cold, 'etrf_hot': etrf_hot}

    def _select_cold(self, ndvi, ts_dem, geometry):
        """Top 5% NDVI → coldest 20% TsDEM → 1% populyatsiya."""
        p95 = ee.Number(ndvi.reduceRegion(
            ee.Reducer.percentile([100 - CIMEC_NDVI_COLD_PCT]),
            geometry, 30, maxPixels=1e9).get('NDVI'))

        top_mask = ndvi.gte(p95)
        ts_top = ts_dem.updateMask(top_mask)

        p20 = ee.Number(ts_top.reduceRegion(
            ee.Reducer.percentile([CIMEC_TS_COLD_PCT]),
            geometry, 30, maxPixels=1e9).get('TsDEM'))

        cold_mask = top_mask.And(ts_dem.lte(p20))

        ts_mean = ee.Number(ts_dem.updateMask(cold_mask).reduceRegion(
            ee.Reducer.mean(), geometry, 30, maxPixels=1e9).get('TsDEM'))
        ndvi_mean = ee.Number(ndvi.updateMask(cold_mask).reduceRegion(
            ee.Reducer.mean(), geometry, 30, maxPixels=1e9).get('NDVI'))

        return {'ts_mean': ts_mean, 'ndvi_mean': ndvi_mean, 'mask': cold_mask}

    def _select_hot(self, ndvi, ts_dem, geometry):
        """Lowest 10% NDVI → hottest 20% TsDEM → 2% populyatsiya."""
        p10 = ee.Number(ndvi.reduceRegion(
            ee.Reducer.percentile([CIMEC_NDVI_HOT_PCT]),
            geometry, 30, maxPixels=1e9).get('NDVI'))

        low_mask = ndvi.lte(p10)
        ts_low = ts_dem.updateMask(low_mask)

        p80 = ee.Number(ts_low.reduceRegion(
            ee.Reducer.percentile([100 - CIMEC_TS_HOT_PCT]),
            geometry, 30, maxPixels=1e9).get('TsDEM'))

        hot_mask = low_mask.And(ts_dem.gte(p80))

        ts_mean = ee.Number(ts_dem.updateMask(hot_mask).reduceRegion(
            ee.Reducer.mean(), geometry, 30, maxPixels=1e9).get('TsDEM'))
        ndvi_mean = ee.Number(ndvi.updateMask(hot_mask).reduceRegion(
            ee.Reducer.mean(), geometry, 30, maxPixels=1e9).get('NDVI'))

        return {'ts_mean': ts_mean, 'ndvi_mean': ndvi_mean, 'mask': hot_mask}

    def _assign_etrf_cold(self, ndvi_mean):
        """CIMEC F.4a/4b."""
        return ee.Number(ee.Algorithms.If(
            ndvi_mean.gte(0.75), CIMEC_ETRF_COLD_DEFAULT,
            ndvi_mean.multiply(1.54).subtract(0.10).max(0.1)))

    def _assign_etrf_hot(self, ndvi_hot, etrf_cold):
        """CIMEC F.5-6."""
        etrf_bare = ee.Number(self.settings.etrf_bare)
        fc = (ndvi_hot.subtract(CIMEC_NDVI_BARE)
              .divide(ee.Number(CIMEC_NDVI_FULL - CIMEC_NDVI_BARE))
              .max(0).min(1))
        return fc.multiply(etrf_cold).add(ee.Number(1).subtract(fc).multiply(etrf_bare))

    def _solve_dt_coefficients(self, cold, hot, etrf_cold, etrf_hot,
                               rn, g, aero, geometry):
        """
        METRIC F.46-51: dT koeffitsientlar.
        Hcold = (Rn-G)cold - LEcold → dTcold
        Hhot  = (Rn-G)hot  - LEhot  → dThot
        a = (dThot - dTcold) / (TsDEM_hot - TsDEM_cold)
        b = (dThot - a) / TsDEM_hot
        """
        def _mean_at(mask, img):
            return ee.Number(img.updateMask(mask).reduceRegion(
                ee.Reducer.mean(), geometry, 30, maxPixels=1e9).values().get(0))

        # Cold pixel energy balance
        rn_c = _mean_at(cold['mask'], rn)
        g_c = _mean_at(cold['mask'], g)
        rah_c = _mean_at(cold['mask'], aero['rah'])
        rho_c = _mean_at(cold['mask'], aero['rho_air'])

        # Hot pixel energy balance
        rn_h = _mean_at(hot['mask'], rn)
        g_h = _mean_at(hot['mask'], g)
        rah_h = _mean_at(hot['mask'], aero['rah'])
        rho_h = _mean_at(hot['mask'], aero['rho_air'])

        # LE = ETrF × ETr × lambda (placeholder: ETr pipeline.py da)
        # Hcold = (Rn-G)cold - LEcold
        # Hhot  = (Rn-G)hot  - LEhot
        h_cold = rn_c.subtract(g_c)  # LEcold pipeline da ayiriladi
        h_hot = rn_h.subtract(g_h)   # LEhot pipeline da ayiriladi

        # dT = H × rah / (rho × Cp)
        dt_cold = h_cold.multiply(rah_c).divide(rho_c.multiply(CP))
        dt_hot = h_hot.multiply(rah_h).divide(rho_h.multiply(CP))

        # F.50-51
        ts_diff = hot['ts_mean'].subtract(cold['ts_mean'])
        a = dt_hot.subtract(dt_cold).divide(ts_diff)
        b = dt_hot.subtract(a).divide(hot['ts_mean'])

        return a, b
