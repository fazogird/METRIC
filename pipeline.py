"""
METRIC-GEE Pipeline — asosiy orkestrator.

Barcha modullarni ketma-ket chaqiradi:
  Inputs → Vegetation → Albedo → Ts → Emissivity → Rn → G →
  Aerodynamics → CIMEC → H(iterative) → ETinst → ETrF → ET24

Tashqaridan faqat: METRICPipeline(cfg).run_monthly(geometry, start, end)
"""
import ee
import math
from datetime import datetime, timedelta
from collections import defaultdict

from metric_gee.config.settings import Settings
from metric_gee.config.constants import *
from metric_gee.inputs.landsat import LandsatLoader
from metric_gee.inputs.era5 import ERA5Loader, prepare_era5_for_pm_daily
from metric_gee.inputs.dem import DEMLoader
from metric_gee.inputs.landuse import LandUseLoader
from metric_gee.core.vegetation import VegetationIndices
from metric_gee.core.albedo import AlbedoCalculator
from metric_gee.core.surface_temp import SurfaceTemperature
from metric_gee.core.emissivity import EmissivityCalculator
from metric_gee.core.net_radiation import NetRadiation
from metric_gee.core.soil_heat_flux import SoilHeatFlux
from metric_gee.core.aerodynamics import Aerodynamics
from metric_gee.core.sensible_heat import SensibleHeatFlux
from metric_gee.core.anchor_selector import AnchorSelector
from metric_gee.core.reference_et import RefETCalculator
from metric_gee.products.instantaneous_et import InstantaneousET
from metric_gee.products.daily_et import DailyET
from metric_gee.utils.solar_geometry import SolarGeometry


class METRICPipeline:
    """METRIC energiya balansi pipeline."""

    def __init__(self, settings=None):
        self.cfg = settings or Settings()
        self._init_modules()

    def _init_modules(self):
        """Barcha modullarni bir marta yaratish."""
        self.landsat = LandsatLoader(self.cfg)
        self.era5 = ERA5Loader(self.cfg)
        self.dem_loader = DEMLoader(self.cfg)
        self.landuse = LandUseLoader(self.cfg)
        self.veg = VegetationIndices(self.cfg)
        self.alb = AlbedoCalculator(self.cfg)
        self.st = SurfaceTemperature(self.cfg)
        self.emiss = EmissivityCalculator()
        self.rn_calc = NetRadiation(self.cfg)
        self.g_calc = SoilHeatFlux(self.cfg)
        self.aero_calc = Aerodynamics(self.cfg)
        self.h_calc = SensibleHeatFlux()
        self.anchor_sel = AnchorSelector(
            cold_lc=self.cfg.cold_lulc_codes or [40],
            hot_lc=self.cfg.hot_lulc_codes   or [60, 50],
            lapse_rate=self.cfg.lapse_rate,
            verbose=True,
        )
        self.et_inst_calc = InstantaneousET()
        self.et_daily = DailyET()
        self.solar_calc = SolarGeometry()
        self.ref_et = RefETCalculator(
            timezone_lon_deg=self.cfg.timezone_lon_deg,
            ref_type='alfalfa')

    # ================================================================
    # PUBLIC API
    # ================================================================

    def run_monthly(self, geometry, start_date, end_date, verbose=True):
        """
        Oylik ET hisoblash — to'liq pipeline.

        Args:
            geometry: ee.Geometry
            start_date: str 'YYYY-MM-DD'
            end_date: str 'YYYY-MM-DD'

        Returns:
            dict: monthly_et, etrf_results, daily_results, diagnostics
        """
        self.geometry = geometry
        self.verbose = verbose

        # DEM — bir marta
        self.dem = self.dem_loader.get_terrain(geometry)
        self.elevation = self.dem.select('elevation')
        self.z_datum = self.dem_loader.get_mean_elevation(geometry)

        # STEP 1: Tasvirlar topish
        images = self._find_images(geometry, start_date, end_date)
        if not images:
            return {'error': 'Tasvir topilmadi'}

        # STEP 2: Har bir tasvir uchun METRIC
        etrf_results = []
        for img_info in images:
            result = self._process_single_image(img_info)
            if result:
                etrf_results.append(result)

        if not etrf_results:
            return {'error': 'Hech bir tasvir ishlamadi'}

        # STEP 3: F.58 oylik ET
        daily_results, monthly_et = self._calc_monthly_et(
            etrf_results, geometry, start_date, end_date)

        return {
            'monthly_et': monthly_et,
            'etrf_results': etrf_results,
            'daily_results': daily_results,
        }

    def run_single(self, geometry, date):
        """Bitta sana uchun METRIC — ETrF xaritasi qaytaradi."""
        self.geometry = geometry
        self.verbose = True
        self.dem = self.dem_loader.get_terrain(geometry)
        self.elevation = self.dem.select('elevation')
        self.z_datum = self.dem_loader.get_mean_elevation(geometry)

        img_info = self._find_best_image(geometry, date)
        if not img_info:
            return None
        return self._process_single_image(img_info)

    # ================================================================
    # STEP 1: TASVIRLAR TOPISH
    # ================================================================

    def _find_images(self, geometry, start, end):
        """Sana oralig'idagi eng yaxshi tasvirlarni topish."""
        l8 = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
            .filterBounds(geometry).filterDate(start, end)
            .filter(ee.Filter.lt('CLOUD_COVER', self.cfg.cloud_cover_max)))
        l9 = (ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
            .filterBounds(geometry).filterDate(start, end)
            .filter(ee.Filter.lt('CLOUD_COVER', self.cfg.cloud_cover_max)))

        # --- Muayyan PATH/ROW bilan cheklash (validatsiya uchun) ---
        if self.cfg.wrs_path is not None and self.cfg.wrs_row is not None:
            path_filter = ee.Filter.And(
                ee.Filter.eq('WRS_PATH', self.cfg.wrs_path),
                ee.Filter.eq('WRS_ROW',  self.cfg.wrs_row)
            )
            l8 = l8.filter(path_filter)
            l9 = l9.filter(path_filter)
            self._log(f"  WRS filter: Path={self.cfg.wrs_path} Row={self.cfg.wrs_row}")

        merged = l8.merge(l9).sort('system:time_start')
        count = merged.size().getInfo()
        self._log(f"Topilgan tasvirlar: {count}")

        if count == 0:
            return []

        image_list = merged.toList(count)

        # Sanalar bo'yicha guruhlash
        date_images = defaultdict(list)
        for i in range(count):
            img = ee.Image(image_list.get(i))
            d = ee.Date(img.get('system:time_start')).format('YYYY-MM-dd').getInfo()
            date_images[d].append((i, img))

        # Har sana uchun eng yaxshi qoplamli tasvirni tanlash
        results = []
        for d, imgs in sorted(date_images.items()):
            best_img = None
            best_count = 0
            for idx, img in imgs:
                cnt = img.select('SR_B4').reduceRegion(
                    ee.Reducer.count(), geometry, 100, maxPixels=1e9
                ).values().get(0)
                cnt_val = ee.Number(cnt).getInfo() or 0
                if cnt_val > best_count:
                    best_count = cnt_val
                    best_img = img
            if best_count > 100:
                results.append({
                    'date': d,
                    'image': best_img,
                    'img_date': ee.Date(best_img.get('system:time_start')),
                    'coverage': best_count,
                })
                self._log(f"  {d}: {best_count} pixels ✅")
            else:
                self._log(f"  {d}: {best_count} pixels ⚠️ skip")

        return results

    def _find_best_image(self, geometry, date):
        """Bitta sana uchun eng yaxshi tasvir."""
        imgs = self._find_images(geometry, date,
            ee.Date(date).advance(16, 'day').format('YYYY-MM-dd').getInfo())
        return imgs[0] if imgs else None

    # ================================================================
    # STEP 2: BITTA TASVIR UCHUN TO'LIQ METRIC
    # ================================================================

    def _process_single_image(self, img_info):
        """Bitta Landsat tasvir uchun to'liq METRIC pipeline."""
        raw_img = img_info['image']
        date_str = img_info['date']
        img_date = img_info['img_date']
        sun_elev = ee.Number(raw_img.get('SUN_ELEVATION'))
        geometry = self.geometry

        self._log(f"\n  --- {date_str} ---")

        try:
            # Mod 0: Input tayyorlash
            image = self.landsat._apply_scale_factors(raw_img)
            image = ee.Image(image.copyProperties(raw_img, raw_img.propertyNames()))
            image = self.landsat._rename_bands(image)
            image = self.landsat._apply_cloud_mask(image)

            # Mod 1: Vegetation
            ndvi, savi, lai = self.veg.calc_all(image)

            # Mod 2: Albedo
            albedo = self.alb.calc_albedo(image)

            # Mod 3: Ts, TsDEM
            ts = self.st.calc_ts(image)
            ts_dem = self.st.calc_ts_dem(ts, self.dem, self.z_datum)

            # Mod 4: Emissivity
            e0, enb = self.emiss.calc_emissivity(lai, ndvi)

            # Mod 5: ERA5 + Rn
            era5_h = self.era5.get_hourly(geometry, img_date)
            solar = self.solar_calc.calc_solar_params(date_str, geometry, sun_elev)
            rn, tau_sw = self.rn_calc.calc_rn(albedo, ts, e0, era5_h, self.dem, solar)

            # Mod 6: G
            g = self.g_calc.calc_g(rn, ts, ndvi, albedo, lai)

            # Mod 7: Aerodynamics
            aero = self.aero_calc.calc_all(era5_h, self.dem, lai, ndvi, albedo, ts)

            # ---- BAND DIAGNOSTIKA ----
            if self.verbose:
                self._print_band_diagnostics(
                    image, ndvi, savi, lai, albedo, ts, ts_dem,
                    e0, rn, g, aero, era5_h, geometry)

            # ETr — Module 5 RS↓ dan (ERA5 ssrd akkumulyatsiya muammosi yo'q)
            rs_down_mjhr = self.rn_calc._calc_rs_down(tau_sw, solar).multiply(0.0036)
            etr_hourly = self._calc_etr_hourly(era5_h, rs_down_mjhr, img_date)

            # Mod 8: CIMEC
            anchor_data, calib = self._run_cimec(
                ndvi, ts_dem, rn, g, aero, etr_hourly, geometry)

            # Mod 9: H — to'liq iteratsiya (dT har qadam yangilanadi)
            H_final = self.h_calc.calc_h_iterative(
                anchor_data, ts_dem, ts, aero, rn, g,
                geometry, verbose=self.verbose)

            # Mod 10: ET
            le, et_inst = self.et_inst_calc.calc(rn, g, H_final, ts)

            # Mod 11: ETrF
            # ── ETr: getInfo() bilan tekshirib olamiz (lazy eval muammosini oldini olish) ──
            etr_dict = etr_hourly.reduceRegion(
                ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
            ).getInfo()
            etr_mean_val = etr_dict.get('ETr', None)

            self._log(f"    [ETr debug] dict keys={list(etr_dict.keys())}  ETr={etr_mean_val}")

            if etr_mean_val is None or etr_mean_val <= 0:
                self._log(f"    ⚠️ ETr={etr_mean_val} — fallback 0.7 mm/hr ishlatiladi")
                etr_mean_val = 0.7   # mm/hr — Idaho July uchun oqilona fallback

            etr_img = ee.Image.constant(etr_mean_val)
            etrf = self.et_daily.calc_etrf(et_inst, etr_img)

            # Statistika — har band alohida (mask propagatsiyasi oldini olish)
            stats = self._get_stats(etrf, et_inst, H_final, rn, g, geometry)
            self._log(f"    a={calib['a_val']:.4f} b={calib['b_val']:.2f}")
            self._log(f"    Rn={stats['Rn']:.0f} G={stats['G']:.0f} H={stats['H']:.0f} W/m²")
            self._log(f"    LE={stats['LE']:.0f} W/m²  ETr={etr_mean_val:.3f} mm/hr")
            self._log(f"    ETrF(tile)={stats['ETrF']:.3f}  ETinst={stats['ETinst']:.3f} mm/hr")

            # ── Qishloq xo'jaligi piksellar uchun ETrF (NDVI ≥ 0.4) ──
            # Sabab: butun tile 95%+ cho'l/tog' — P50=0.
            # Sug'oriladigan ekinlar NDVI≥0.4 bo'ladi (o'simlik ko'p).
            agri_ndvi_thresh = 0.4
            agri_mask = ndvi.gte(agri_ndvi_thresh)
            etrf_agri = etrf.updateMask(agri_mask)

            agri_s = etrf_agri.reduceRegion(
                ee.Reducer.mean()
                          .combine(ee.Reducer.percentile([25, 50, 75, 90]),
                                   sharedInputs=True),
                geometry, 200, maxPixels=1e9
            ).getInfo()

            etrf_agri_mean = agri_s.get('ETrF_mean') or agri_s.get('ETrF') or 0.0
            etrf_agri_p50  = agri_s.get('ETrF_p50')  or etrf_agri_mean
            etrf_agri_p75  = agri_s.get('ETrF_p75')  or etrf_agri_p50

            self._log(
                f"    ETrF agri(NDVI≥{agri_ndvi_thresh}): "
                f"mean={etrf_agri_mean:.3f}  "
                f"P25={agri_s.get('ETrF_p25',0):.3f}  "
                f"P50={etrf_agri_p50:.3f}  "
                f"P75={etrf_agri_p75:.3f}  "
                f"P90={agri_s.get('ETrF_p90',0):.3f}"
            )

            # Butun tile uchun ham percentile (tashxis uchun)
            tile_pct = etrf.reduceRegion(
                ee.Reducer.percentile([75, 90, 95]),
                geometry, 500, maxPixels=1e9
            ).getInfo()
            self._log(
                f"    ETrF tile: "
                f"P75={tile_pct.get('ETrF_p75',0):.3f}  "
                f"P90={tile_pct.get('ETrF_p90',0):.3f}  "
                f"P95={tile_pct.get('ETrF_p95',0):.3f}"
            )

            # etrf_mean = qishloq xo'jaligi piksellar o'rtachasi (oylik ET uchun)
            return {
                'date': date_str,
                'etrf_mean':      etrf_agri_mean,   # ← NDVI≥0.4 o'rtacha
                'etrf_agri_p75':  etrf_agri_p75,
                'etrf_tile_p90':  tile_pct.get('ETrF_p90', 0),
                'etrf_image':     etrf,
                'stats':          stats,
                'etr_hourly_val': etr_mean_val,
            }

        except Exception as e:
            self._log(f"    ⚠️ Xatolik: {e}")
            return None

    # ================================================================
    # ETr HISOBLASH — RefETCalculator orqali
    # ================================================================

    def _calc_etr_hourly(self, era5_h, rs_down_mjhr, img_date):
        """
        Soatlik ETr — RefETCalculator bilan.
        Rs = Module 5 dan (ERA5 ssrd akkumulyatsiya muammosi yo'q).
        """
        T_air = era5_h.select('temperature_2m').subtract(273.15).rename('T_air')
        P_kPa = era5_h.select('sp_kpa').rename('P_kPa')
        u2 = era5_h.select('wind_speed_10m').multiply(0.748).rename('u2')
        ea = era5_h.select('ea_kpa').rename('ea')

        era5_pm = (T_air.addBands(P_kPa).addBands(u2).addBands(ea)
            .addBands(rs_down_mjhr.rename('Rs'))
            .set('system:time_start', img_date.millis()))

        result = self.ref_et.calculate(era5_pm, self.elevation, mode='hourly')
        return result.select('ETr')

    # ================================================================
    # CIMEC KALIBRLASH
    # ================================================================

    def _run_cimec(self, ndvi, ts_dem, rn, g, aero, etr_hourly, geometry):
        """
        CIMEC — LULC asosida cold/hot piksel tanlash + energiya balansi.

        1) AnchorSelector → LULC cascade bilan cold/hot maskalar
        2) Anchor nuqtada Rn, G, rah, rho o'rtachasini olish
        3) LE anchor → H anchor hisoblash
        """
        # ── 1) LULC asosida anchor tanlash ──────────────────────
        sel = self.anchor_sel.select(ndvi, ts_dem, geometry)

        cold_mask = sel['cold_mask']
        hot_mask  = sel['hot_mask']
        cold_ts   = sel['cold_ts']
        hot_ts    = sel['hot_ts']

        self._log(f"  Anchor plan: {sel['plan']}")
        self._log(f"  Cold: Ts={cold_ts.getInfo():.2f} K  "
                  f"NDVI={sel['cold_ndvi'].getInfo():.3f}")
        self._log(f"  Hot:  Ts={hot_ts.getInfo():.2f} K  "
                  f"NDVI={sel['hot_ndvi'].getInfo():.3f}")

        # ── 2) Anchor nuqtalarda energiya parametrlari ───────────
        def anchor_mean(img, mask, band):
            res = img.updateMask(mask).reduceRegion(
                ee.Reducer.mean(), geometry, 100, maxPixels=1e9
            ).getInfo()
            for key in [band, f'{band}_mean']:
                if key in res and res[key] is not None:
                    return ee.Number(res[key])
            vals = [v for v in res.values() if v is not None]
            return ee.Number(vals[0] if vals else 0.0)

        cold_rn  = anchor_mean(rn,             cold_mask, 'Rn')
        cold_g   = anchor_mean(g,              cold_mask, 'G')
        hot_rn   = anchor_mean(rn,             hot_mask,  'Rn')
        hot_g    = anchor_mean(g,              hot_mask,  'G')
        cold_rah = anchor_mean(aero['rah'],    cold_mask, 'rah')
        hot_rah  = anchor_mean(aero['rah'],    hot_mask,  'rah')
        cold_rho = anchor_mean(aero['rho_air'],cold_mask, 'rho_air')
        hot_rho  = anchor_mean(aero['rho_air'],hot_mask,  'rho_air')

        self._log(f"  Cold anchor: Rn={cold_rn.getInfo():.1f}  "
                  f"G={cold_g.getInfo():.1f}  "
                  f"rah={cold_rah.getInfo():.1f}  "
                  f"rho={cold_rho.getInfo():.4f}")
        self._log(f"  Hot  anchor: Rn={hot_rn.getInfo():.1f}  "
                  f"G={hot_g.getInfo():.1f}  "
                  f"rah={hot_rah.getInfo():.1f}  "
                  f"rho={hot_rho.getInfo():.4f}")

        # ── 3) ETr scalar ────────────────────────────────────────
        etr_val = ee.Number(etr_hourly.reduceRegion(
            ee.Reducer.mean(), geometry, 5000, maxPixels=1e9).get('ETr'))

        # ── 4) LE anchor (F.53) ──────────────────────────────────
        # LEcold = ETrF_cold × ETr × λ / 3600  (ETrF_cold = 1.05)
        lam_cold = cold_ts.subtract(273.15).multiply(-0.00236).add(2.501).multiply(1e6)
        le_cold  = ee.Number(CIMEC_ETRF_COLD_DEFAULT).multiply(etr_val) \
                             .multiply(lam_cold).divide(3600)

        # LEhot = ETrF_bare × ETr × λ / 3600  (ETrF_bare = 0.05)
        lam_hot = hot_ts.subtract(273.15).multiply(-0.00236).add(2.501).multiply(1e6)
        le_hot  = ee.Number(self.cfg.etrf_bare).multiply(etr_val) \
                            .multiply(lam_hot).divide(3600)

        # ── 5) H anchor = (Rn - G) - LE ─────────────────────────
        h_cold = cold_rn.subtract(cold_g).subtract(le_cold)
        h_hot  = hot_rn.subtract(hot_g).subtract(le_hot)

        self._log(f"  LE_cold={le_cold.getInfo():.1f}  H_cold={h_cold.getInfo():.1f} W/m²")
        self._log(f"  LE_hot={le_hot.getInfo():.1f}   H_hot={h_hot.getInfo():.1f} W/m²")

        # ── 6) Boshlang'ich a, b (log uchun) ─────────────────────
        dt_c = h_cold.multiply(cold_rah).divide(cold_rho.multiply(CP))
        dt_h = h_hot.multiply(hot_rah).divide(hot_rho.multiply(CP))
        ts_d = hot_ts.subtract(cold_ts)
        a_init = dt_h.subtract(dt_c).divide(ts_d)
        b_init = dt_h.subtract(a_init.multiply(hot_ts))

        calib = {
            'a_val': a_init.getInfo(),
            'b_val': b_init.getInfo(),
        }

        anchor_data = {
            'h_cold':    h_cold,
            'h_hot':     h_hot,
            'cold_mask': cold_mask,
            'hot_mask':  hot_mask,
            'cold_ts':   cold_ts,
            'hot_ts':    hot_ts,
            'cold_rah':  cold_rah,
            'hot_rah':   hot_rah,
            'cold_rho':  cold_rho,
            'hot_rho':   hot_rho,
        }

        return anchor_data, calib

    # ================================================================
    # H ITERATSIYA
    # ================================================================

    def _run_h_iteration(self, a, b, ts_dem, ts, aero, rn, g):
        """
        Monin-Obukhov iterativ H hisoblash (10 iteratsiya).

        Fizik cheklovlar (Allen 2007 + ASCE-EWRI):
          rah >= 5 s/m  (aerodinamik razshilik minimum)
          |L| >= 2 m    (Obukhov uzunligi singularyatdan saqlash)
          H  <= (Rn-G)  (LE manfiy bo'la olmaydi, fizik limit)
          H  >= -0.5*(Rn-G) (kechqurun kondensatsiya uchun kichik marja)
        """
        dT    = ts_dem.multiply(a).add(b).rename('dT')
        u200  = aero['u200']
        zom   = aero['zom']
        u_star = aero['u_star']
        rah   = aero['rah']
        rho_air = aero['rho_air']

        # Fizik chegaralar konstantasi
        rah_min  = ee.Image.constant(5.0)    # s/m — Allen: rah ≥ 5
        rah_max  = ee.Image.constant(500.0)  # s/m
        L_min_abs = 2.0                       # |L| ≥ 2 m

        for i in range(10):
            H = rho_air.multiply(CP).multiply(dT).divide(rah)

            # H=0 bo'lganda L→∞ → stable limit
            H_safe = H.where(H.abs().lt(1.0), ee.Image.constant(1.0))

            L = (rho_air.multiply(CP).multiply(u_star.pow(3)).multiply(ts)
                 .divide(ee.Image.constant(VON_KARMAN * GRAVITY).multiply(H_safe))
                 .multiply(-1))

            # |L| >= L_min_abs, belgisi saqlansin
            is_stable  = L.gt(0)
            L = L.abs().max(L_min_abs).multiply(
                is_stable.multiply(2).subtract(1).multiply(-1)  # stable→+, unstable→-
            )
            # Mantiqan to'g'rilash: unstable L < 0, stable L > 0
            L = L.where(is_stable, L.abs())
            L = L.where(is_stable.Not(), L.abs().multiply(-1))
            L = L.clamp(-2000, 2000)

            # Beqarorlik/barqarorlik korreksiya funksiyalari (Paulson 1970)
            x200 = ee.Image.constant(1).subtract(ee.Image.constant(16*Z_BLEND).divide(L)).pow(0.25)
            x2   = ee.Image.constant(1).subtract(ee.Image.constant(16*Z2).divide(L)).pow(0.25)
            x01  = ee.Image.constant(1).subtract(ee.Image.constant(16*Z1).divide(L)).pow(0.25)

            psi_m   = (x200.add(1).divide(2).log().multiply(2)
                       .add(x200.pow(2).add(1).divide(2).log())
                       .subtract(x200.atan().multiply(2)).add(math.pi/2))
            psi_h2  = x2.pow(2).add(1).divide(2).log().multiply(2)
            psi_h01 = x01.pow(2).add(1).divide(2).log().multiply(2)

            # Barqaror: Webb (1970) linear korreksiya
            psi_m_s   = ee.Image.constant(-5 * Z_BLEND).divide(L)
            psi_h2_s  = ee.Image.constant(-5 * Z2).divide(L)
            psi_h01_s = ee.Image.constant(-5 * Z1).divide(L)

            psi_m   = psi_m.where(is_stable, psi_m_s)
            psi_h2  = psi_h2.where(is_stable, psi_h2_s)
            psi_h01 = psi_h01.where(is_stable, psi_h01_s)

            # Psi cheklov: psi_m va psi_h juda katta bo'lmasin
            psi_m   = psi_m.clamp(-10, 10)
            psi_h2  = psi_h2.clamp(-10, 10)
            psi_h01 = psi_h01.clamp(-10, 10)

            denom_m = (ee.Image.constant(Z_BLEND).divide(zom).log()
                       .subtract(psi_m).max(0.5))  # manfiy denominatorsdan saqlash
            u_star = (u200.multiply(VON_KARMAN).divide(denom_m).max(0.01))

            rah = (ee.Image.constant(math.log(Z2/Z1))
                   .subtract(psi_h2).add(psi_h01)
                   .divide(u_star.multiply(VON_KARMAN)))
            rah = rah.where(rah.lt(0), rah_min)   # manfiy rah → minimum
            rah = rah.max(rah_min).min(rah_max)    # [5, 500] s/m

        H_final = rho_air.multiply(CP).multiply(dT).divide(rah).rename('H')

        # LE manfiy bo'lmasin (Rn-G fizik limiti)
        rn_g = rn.subtract(g)
        H_final = H_final.min(rn_g)                  # H ≤ Rn-G → LE ≥ 0
        H_final = H_final.max(rn_g.multiply(-0.3))   # H ≥ -0.3*(Rn-G) (kechki kondensatsiya)

        return H_final

    # ================================================================
    # STEP 3: F.58 OYLIK ET
    # ================================================================

    def _calc_monthly_et(self, etrf_results, geometry, start, end):
        """F.58: ETperiod = Σ(ETrFi × ETr24i)"""
        start_dt = datetime.strptime(start, '%Y-%m-%d')
        end_dt = datetime.strptime(end, '%Y-%m-%d')
        n_days = (end_dt - start_dt).days + 1

        # Tasvir sanalari va ETrF
        img_dates = [datetime.strptime(r['date'], '%Y-%m-%d') for r in etrf_results]
        img_etrfs = [r['etrf_mean'] for r in etrf_results]

        self._log(f"\n  ETrF nuqtalar: {len(img_dates)} ta")
        for d, e in zip(img_dates, img_etrfs):
            self._log(f"    {d.strftime('%b %d')} → ETrF = {e:.3f}")

        # Kunlik hisob
        daily_results = []
        et_total = 0
        etr_total = 0

        for day in range(n_days):
            current = start_dt + timedelta(days=day)
            date_str = current.strftime('%Y-%m-%d')

            # ETrF interpolatsiya
            etrf_today = self._interpolate_etrf(current, img_dates, img_etrfs)

            # Kunlik ETr — RefETCalculator
            try:
                era5_day = prepare_era5_for_pm_daily(geometry, date_str)
                result = self.ref_et.calculate(era5_day, self.elevation, mode='daily')
                etr_today = ee.Number(result.select('ETr').reduceRegion(
                    ee.Reducer.mean(), geometry, 1000, maxPixels=1e9
                ).get('ETr')).getInfo()
                etr_today = max(etr_today or 0.1, 0.1)
            except:
                etr_today = 7.0

            et_today = etrf_today * etr_today
            et_total += et_today
            etr_total += etr_today

            is_img = "◄" if current in img_dates else ""
            daily_results.append({
                'date': date_str, 'etrf': etrf_today,
                'etr': etr_today, 'et': et_today})

            self._log(f"  {date_str}  {etrf_today:6.3f}  {etr_today:7.2f}  {et_today:7.2f}  {is_img}")

        self._log(f"  SUMMA: ETr={etr_total:.1f} ET={et_total:.1f} mm/oy")

        return daily_results, {
            'et_monthly': et_total,
            'etr_monthly': etr_total,
            'etrf_weighted': et_total / max(etr_total, 0.1),
        }

    # ================================================================
    # YORDAMCHI FUNKSIYALAR
    # ================================================================

    @staticmethod
    def _interpolate_etrf(target, dates, etrfs):
        """Linear interpolatsiya — flat extrapolation chegaralarda."""
        if len(dates) == 1:
            return etrfs[0]
        if target <= dates[0]:
            return etrfs[0]
        if target >= dates[-1]:
            return etrfs[-1]
        for i in range(len(dates)-1):
            if dates[i] <= target <= dates[i+1]:
                total = (dates[i+1]-dates[i]).days
                frac = (target-dates[i]).days / max(total, 1)
                return etrfs[i] + frac * (etrfs[i+1]-etrfs[i])
        return etrfs[-1]

    @staticmethod
    def _percentile(img, band, pct, geometry):
        """
        ee.Image dan percentile olish.
        GEE versiyasiga qarab kalit 'BAND_pN' yoki 'BAND' shaklida bo'lishi mumkin.
        """
        result = img.reduceRegion(
            ee.Reducer.percentile([pct]), geometry, 30, maxPixels=1e9).getInfo()

        if not result:
            raise ValueError(f"_percentile: '{band}_p{pct}' — bo'sh natija (barcha piksellar maskalangan?)")

        # GEE ba'zan '_pN' qo'shadi, ba'zan qo'shmaydi
        for key in [f'{band}_p{pct}', band]:
            if key in result and result[key] is not None:
                return ee.Number(result[key])

        raise ValueError(
            f"_percentile: '{band}_p{pct}' yoki '{band}' kaliti topilmadi. "
            f"Mavjud kalitlar: {list(result.keys())}"
        )

    @staticmethod
    def _masked_mean(img, mask, band, geometry):
        """Maskalangan piksellardan o'rtacha. Kalit yo'q bo'lsa None."""
        return ee.Number(img.updateMask(mask).reduceRegion(
            ee.Reducer.mean(), geometry, 30, maxPixels=1e9).get(band, None))

    def _get_stats(self, etrf, et_inst, h, rn, g, geometry):
        """
        O'rtacha statistika — har band ALOHIDA reduce qilinadi.
        Sabab: addBands() mask intersection qiladi va bitta
        maskalangan band barcha bandlarni nolga tushirishi mumkin.
        """
        scale = 100   # 100m — tez hisoblash uchun (30m juda sekin)

        def _mean(img, band):
            try:
                res = img.reduceRegion(
                    ee.Reducer.mean(), geometry, scale, maxPixels=1e9
                ).getInfo()
                # Band nomini avto-topish (ba'zan GEE band nomi o'zgarishi mumkin)
                if band in res:
                    v = res[band]
                else:
                    v = list(res.values())[0] if res else None
                return v if v is not None else 0.0
            except Exception as ex:
                return 0.0

        # LE = Rn - G - H (alohida hisoblash)
        le_img = rn.subtract(g).subtract(h).rename('LE')

        return {
            'Rn':     _mean(rn,     'Rn'),
            'G':      _mean(g,      'G'),
            'H':      _mean(h,      'H'),
            'LE':     _mean(le_img, 'LE'),
            'ETinst': _mean(et_inst,'ETinst'),
            'ETrF':   _mean(etrf,   'ETrF'),
        }

    # ================================================================
    # BAND DIAGNOSTIKA
    # ================================================================

    def _print_band_diagnostics(self, image, ndvi, savi, lai, albedo,
                                 ts, ts_dem, e0, rn, g, aero, era5_h, geometry):
        """
        Barcha muhim band/parametrlarning min/mean/max ni chop etish.
        Muammoni topish uchun — fizik chegaralarni tekshiring.
        """
        scale = 100

        def mmm(img, bands, label):
            """Min / Mean / Max — chiroyli ko'rinish."""
            combined = img.select(bands[0])
            for b in bands[1:]:
                combined = combined.addBands(img.select(b))
            s = combined.reduceRegion(
                ee.Reducer.min().combine(ee.Reducer.mean(), sharedInputs=True)
                          .combine(ee.Reducer.max(), sharedInputs=True),
                geometry, scale, maxPixels=1e9
            ).getInfo()

            self._log(f"    [{label}]")
            for b in bands:
                mn  = s.get(f'{b}_min',  s.get(f'{b}', '?'))
                avg = s.get(f'{b}_mean', '?')
                mx  = s.get(f'{b}_max',  '?')
                if all(isinstance(v, (int, float)) for v in [mn, avg, mx]):
                    self._log(f"      {b:<20s}  min={mn:9.4f}  mean={avg:9.4f}  max={mx:9.4f}")
                else:
                    self._log(f"      {b:<20s}  {s}")

        self._log("  ┌─ BAND DIAGNOSTIKA ─────────────────────────────────────────")

        # --- Landsat SR ---
        try:
            mmm(image, ['B_BLUE','B_GREEN','B_RED','B_NIR','B_SWIR1','B_SWIR2'],
                "Landsat SR (surface reflectance, 0-1)")
            mmm(image, ['B_THERMAL'], "Landsat Ts raw (K, birinchi qadam)")
        except Exception as e:
            self._log(f"      Landsat SR: xato — {e}")

        # --- Vegetation ---
        try:
            mmm(ndvi,   ['NDVI'],  "Vegetation: NDVI (-1..1)")
            mmm(savi,   ['SAVI'],  "Vegetation: SAVI")
            mmm(lai,    ['LAI'],   "Vegetation: LAI (m²/m²)")
        except Exception as e:
            self._log(f"      Vegetation: xato — {e}")

        # --- Albedo ---
        try:
            mmm(albedo, ['albedo'], "Albedo (0..0.4)")
        except Exception as e:
            self._log(f"      Albedo: xato — {e}")

        # --- Surface Temperature ---
        try:
            mmm(ts,     ['Ts'],    "Ts — actual (K)")
            mmm(ts_dem, ['TsDEM'], "TsDEM — lapse corrected (K)")
            mmm(e0,     ['e0'],    "Emissivity e0 (0.95..0.99)")
        except Exception as e:
            self._log(f"      Ts/Emiss: xato — {e}")

        # --- Energy Balance ---
        try:
            mmm(rn, ['Rn'], "Rn — net radiation (W/m²)")
            mmm(g,  ['G'],  "G  — soil heat flux (W/m²)")
        except Exception as e:
            self._log(f"      Rn/G: xato — {e}")

        # --- Aerodynamics ---
        try:
            mmm(aero['u200'],    ['u200'],    "u200 — wind 200m (m/s)")
            mmm(aero['zom'],     ['zom'],     "zom  — roughness (m)")
            mmm(aero['u_star'],  ['u_star'],  "u*   — friction vel (m/s)")
            mmm(aero['rah'],     ['rah'],     "rah  — aero resist (s/m) ← KEY!")
            mmm(aero['rho_air'], ['rho_air'], "rho  — air density (kg/m³)")
        except Exception as e:
            self._log(f"      Aero: xato — {e}")

        # --- ERA5 (scalar, o'rtacha qiymat) ---
        try:
            era5_vars = [
                'temperature_2m', 'dewpoint_temperature_2m',
                'wind_speed_10m', 'sp_kpa', 'ea_kpa',
            ]
            s = era5_h.select(era5_vars).reduceRegion(
                ee.Reducer.mean(), geometry, 10000, maxPixels=1e9
            ).getInfo()
            self._log("    [ERA5 hourly — o'rtacha]")
            labels = {
                'temperature_2m':         'T2m (K)',
                'dewpoint_temperature_2m':'Td  (K)',
                'wind_speed_10m':         'WS10 (m/s)',
                'sp_kpa':                 'P (kPa)',
                'ea_kpa':                 'ea (kPa)',
            }
            for k, lbl in labels.items():
                v = s.get(k, '?')
                if isinstance(v, (int, float)):
                    self._log(f"      {lbl:<22s}  {v:.4f}")
        except Exception as e:
            self._log(f"      ERA5: xato — {e}")

        self._log("  └────────────────────────────────────────────────────────────")

    def _log(self, msg):
        if self.verbose:
            print(msg)
