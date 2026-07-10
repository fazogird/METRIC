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
                # LULC cropland ustidagi piksellik bulut tekshiruvi
                if not self._cloud_precheck(best_img, geometry):
                    self._log(f"  {d}: {best_count} pixels — cloud precheck fail → skip")
                    continue
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
            # F.7/F.8: DEM berilsa slope/aspect tuzatma (ee.Image),
            #          berilmasa tekis yer (ee.Number)
            # solar = self.solar_calc.calc_solar_params(
            #     date_str, geometry, sun_elev, dem=self.dem)
            
            solar = self.solar_calc.calc_solar_params(
                date_str, geometry, sun_elev,
                dem=self.dem,
                img_date=img_date,
                timezone_lon_deg=self.cfg.timezone_lon_deg)
            
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

            # Mod 8: CIMEC — F.7 (albedo), F.8 (Tfac), F.4a, F.5/F.6
            anchor_data, calib = self._run_cimec(
                ndvi, ts_dem, rn, g, aero, etr_hourly, geometry,
                albedo   = albedo,
                solar    = solar,
                date_str = date_str)

            # Mod 9: H — to'liq iteratsiya (dT har qadam yangilanadi)
            H_final = self.h_calc.calc_h_iterative(
                anchor_data, ts_dem, ts, aero, rn, g,
                geometry, verbose=self.verbose)

            # Mod 10: ET
            le, et_inst = self.et_inst_calc.calc(rn, g, H_final, ts)

            # Mod 11: ETrF  (faqat alfalfa ETr bilan)
            etrf = self.et_daily.calc_etrf(et_inst, etr_hourly.select('ETr'))

            # Statistika — har band alohida (mask propagatsiyasi oldini olish)
            stats = self._get_stats(etrf, et_inst, H_final, rn, g, geometry)
            ref_vals = etr_hourly.reduceRegion(
                ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
            ).getInfo()
            etr_log = ref_vals.get('ETr', 0) or 0
            eto_log = ref_vals.get('ETo', 0) or 0
            self._log(f"    a={calib['a_val']:.4f} b={calib['b_val']:.2f}")
            self._log(f"    Rn={stats['Rn']:.0f} G={stats['G']:.0f} H={stats['H']:.0f} W/m²")
            self._log(f"    LE={stats['LE']:.0f} W/m²  ETr≈{etr_log:.3f} ETo≈{eto_log:.3f} mm/hr")


            # ── Butun tile ETrF — _get_stats da allaqachon hisoblangan ──
            etrf_tile_mean = stats.get('ETrF', 0.0)
            self._log(f"    ETrF(tile mean)={etrf_tile_mean:.3f}")

            # Mod 11.5: Crad — F.56 qiyalik tuzatma koeffitsienti
            crad = self.et_daily.calc_crad(self.dem, solar)
            crad_mean = crad.reduceRegion(
                ee.Reducer.mean(), geometry, 300, maxPixels=1e9
            ).getInfo().get('Crad', 1.0) or 1.0
            self._log(f"    Crad(tile mean)={crad_mean:.3f}")

            return {
                'date':           date_str,
                'etrf_mean':      etrf_tile_mean,   # ← butun tile, NDVI mask yo'q
                'etrf_image':     etrf,
                'crad_image':     crad,             # ← F.56 qiyalik tuzatma (ee.Image)
                'crad_mean':      crad_mean,
                'stats':          stats,
                'etr_hourly_val': etr_log,
                'eto_hourly_val': eto_log,
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
        return result.select(['ETr', 'ETo'])

    # ================================================================
    # CIMEC KALIBRLASH
    # ================================================================

    def _run_cimec(self, ndvi, ts_dem, rn, g, aero, etr_hourly, geometry,
                   albedo=None, solar=None, date_str=None):
        """
        CIMEC — Allen et al. (2013) to'liq implementatsiya.

        1) Precip_60, ETr_60 olish  → Tfac [F.8]
        2) AnchorSelector → LULC cascade + ±0.2K + albedo [F.7] + Tfac [F.8]
        3) F.4a → ETrFcold,  F.5/F.6 → ETrFhot
        4) LE anchor → H anchor → a, b regression
        """
        # ── 1) F.8 uchun P₆₀, ETr₆₀ ────────────────────────────
        precip_60_val = None
        etr_60_val    = None
        if date_str is not None:
            try:
                p_img = self.era5.get_precip_sum(geometry, date_str, 60)
                precip_60_val = float(
                    p_img.reduceRegion(
                        ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
                    ).getInfo().get('precip_sum_mm') or 0.0)

                e_img = self.era5.get_etr_sum(geometry, date_str, 60)
                etr_60_val = float(
                    e_img.reduceRegion(
                        ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
                    ).getInfo().get('etr_sum_mm') or 1.0)

                self._log(f"  P₆₀={precip_60_val:.1f} mm  ETr₆₀={etr_60_val:.1f} mm")
            except Exception as ex:
                self._log(f"  ⚠️ P₆₀/ETr₆₀ xato: {ex}")

        # ── 2) LULC asosida anchor tanlash (F.7, F.8 bilan) ─────
        sun_elev_val = None
        if solar is not None:
            try:
                sun_elev_val = float(solar['sun_elevation'].getInfo())
            except Exception:
                pass

        sel = self.anchor_sel.select(
            ndvi, ts_dem, geometry,
            albedo        = albedo,
            sun_elev      = sun_elev_val,
            precip_60_val = precip_60_val,
            etr_60_val    = etr_60_val,
        )

        # Mask EPSG:4326 (ESA WorldCover) → Landsat UTM (scale=300m).
        # scale=30 → memory error; scale=300 → 1.2M piksel, muammo yo'q.
        landsat_proj     = rn.projection()
        cold_mask        = sel['cold_mask'].reproject(crs=landsat_proj, scale=300)
        hot_mask         = sel['hot_mask'].reproject(crs=landsat_proj, scale=300)
        cold_ts          = sel['cold_ts']
        hot_ts           = sel['hot_ts']                # fizik (λ uchun)
        hot_ts_corrected = sel['hot_ts_corrected']      # F.8 tuzatilgan (regression uchun)
        cold_ndvi_val    = float(sel['cold_ndvi'].getInfo())
        hot_ndvi_val     = float(sel['hot_ndvi'].getInfo())
        tfac_val         = float(sel['tfac'].getInfo())

        self._log(f"  Anchor plan: {sel['plan']}")
        self._log(f"  Cold: Ts={cold_ts.getInfo():.2f} K  NDVI={cold_ndvi_val:.3f}")
        self._log(f"  Hot:  Ts_phys={hot_ts.getInfo():.2f} K  "
                  f"Tfac={tfac_val:.2f} K  "
                  f"Ts_corr={hot_ts_corrected.getInfo():.2f} K  "
                  f"NDVI={hot_ndvi_val:.3f}")

        # ── 3) F.4a: ETrFcold → F.5/F.6: ETrFhot ───────────────
        # F.4a: ETrFcold = 1.54·NDVI₁ - 0.10  (NDVI₁ < 0.75)
        #        ETrFcold = 1.05              (NDVI₁ ≥ 0.75)
        if cold_ndvi_val >= 0.75:
            etrf_cold = CIMEC_ETRF_COLD_DEFAULT          # 1.05  (F.4b)
            self._log(f"  F.4b ETrFcold = 1.05  (NDVI₁={cold_ndvi_val:.3f} ≥ 0.75)")
        else:
            etrf_cold = max(0.1, 1.54 * cold_ndvi_val - 0.10)   # F.4a
            self._log(f"  F.4a ETrFcold = 1.54×{cold_ndvi_val:.3f} - 0.10 = {etrf_cold:.3f}")

        # F.5: fc = (NDVIhot - NDVIbare) / (NDVIfull - NDVIbare)
        fc = max(0.0, min(1.0,
            (hot_ndvi_val - CIMEC_NDVI_BARE) / (CIMEC_NDVI_FULL - CIMEC_NDVI_BARE)
        ))
        # F.6: ETrFhot = fc · ETrFcold + (1 - fc) · ETrFbare
        etrf_hot = fc * etrf_cold + (1.0 - fc) * self.cfg.etrf_bare
        self._log(f"  F.5/F.6 fc={fc:.3f}  ETrFhot={etrf_hot:.3f}  "
                  f"(NDVIhot={hot_ndvi_val:.3f}  ETrFbare={self.cfg.etrf_bare})")

        # ── 2) Anchor nuqtalarda energiya parametrlari ───────────
        def anchor_mean(img, mask, band):
            # Mask endi Landsat UTM 300m da (yuqorida reproject qilingan).
            for scale in [300, 1000]:
                res = img.updateMask(mask).reduceRegion(
                    ee.Reducer.mean(), geometry, scale, maxPixels=1e9
                ).getInfo()
                for key in [band, f'{band}_mean']:
                    if key in res and res[key] is not None:
                        return ee.Number(res[key])
            # Last resort: maska holda butun qism
            res = img.reduceRegion(
                ee.Reducer.mean(), geometry, 300, maxPixels=1e9
            ).getInfo()
            for key in [band, f'{band}_mean']:
                if key in res and res[key] is not None:
                    return ee.Number(res[key])
            return ee.Number(0.0)

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
        _etr_info = etr_hourly.reduceRegion(
            ee.Reducer.mean(), geometry, 5000, maxPixels=1e9).getInfo()
        _etr_raw = (_etr_info or {}).get('ETr')
        if _etr_raw is None:
            self._log("  ⚠️ ETr null (barcha piksellar maskalangan?) — 0.001 ishlatildi")
        etr_val = ee.Number(_etr_raw if _etr_raw is not None else 0.001)

        # ── 4) LE anchor (F.53) ──────────────────────────────────
        # F.4a/b: ETrFcold  (yuqorida hisoblangan Python float)
        lam_cold = cold_ts.subtract(273.15).multiply(-0.00236).add(2.501).multiply(1e6)
        le_cold  = ee.Number(etrf_cold).multiply(etr_val).multiply(lam_cold).divide(3600)

        # F.5/F.6: ETrFhot  (yuqorida hisoblangan Python float)
        # Fizik hot_ts ishlatiladi — λ = f(T_fizik)
        lam_hot = hot_ts.subtract(273.15).multiply(-0.00236).add(2.501).multiply(1e6)
        le_hot  = ee.Number(etrf_hot).multiply(etr_val).multiply(lam_hot).divide(3600)

        # ── 5) H anchor = (Rn - G) - LE ─────────────────────────
        h_cold = cold_rn.subtract(cold_g).subtract(le_cold)
        h_hot  = hot_rn.subtract(hot_g).subtract(le_hot)

        self._log(f"  LE_cold={le_cold.getInfo():.1f}  H_cold={h_cold.getInfo():.1f} W/m²")
        self._log(f"  LE_hot={le_hot.getInfo():.1f}   H_hot={h_hot.getInfo():.1f} W/m²")

        # ── 6) a, b regression  (F.8: hot_ts_corrected ishlatiladi) ──
        dt_c   = h_cold.multiply(cold_rah).divide(cold_rho.multiply(CP))
        dt_h   = h_hot.multiply(hot_rah).divide(hot_rho.multiply(CP))
        ts_d   = hot_ts_corrected.subtract(cold_ts)     # ← F.8 tuzatilgan ΔTs
        a_init = dt_h.subtract(dt_c).divide(ts_d)
        b_init = dt_h.subtract(a_init.multiply(hot_ts_corrected))  # ← F.8

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
            'hot_ts':    hot_ts_corrected,   # F.8: regression uchun tuzatilgan
            'cold_rah':  cold_rah,
            'hot_rah':   hot_rah,
            'cold_rho':  cold_rho,
            'hot_rho':   hot_rho,
        }

        return anchor_data, calib

 
    # ================================================================
    # STEP 3: F.58 OYLIK ET
    # ================================================================

    # def _calc_monthly_et(self, etrf_results, geometry, start, end):
    #     """F.58: ETperiod = Σ(ETrFi × ETr24i)"""
    #     start_dt = datetime.strptime(start, '%Y-%m-%d')
    #     end_dt = datetime.strptime(end, '%Y-%m-%d')
    #     n_days = (end_dt - start_dt).days + 1

    #     # Tasvir sanalari va ETrF
    #     img_dates = [datetime.strptime(r['date'], '%Y-%m-%d') for r in etrf_results]
    #     img_etrfs = [r['etrf_mean'] for r in etrf_results]

    #     self._log(f"\n  ETrF nuqtalar: {len(img_dates)} ta")
    #     for d, e in zip(img_dates, img_etrfs):
    #         self._log(f"    {d.strftime('%b %d')} → ETrF = {e:.3f}")

    #     # Kunlik hisob
    #     daily_results = []
    #     et_total = 0
    #     etr_total = 0

    #     for day in range(n_days):
    #         current = start_dt + timedelta(days=day)
    #         date_str = current.strftime('%Y-%m-%d')

    #         # ETrF interpolatsiya
    #         etrf_today = self._interpolate_etrf(current, img_dates, img_etrfs)

    #         # Kunlik ETr — RefETCalculator
    #         try:
    #             era5_day = prepare_era5_for_pm_daily(geometry, date_str)
    #             result = self.ref_et.calculate(era5_day, self.elevation, mode='daily')
    #             etr_today = ee.Number(result.select('ETr').reduceRegion(
    #                 ee.Reducer.mean(), geometry, 1000, maxPixels=1e9
    #             ).get('ETr')).getInfo()
    #             etr_today = max(etr_today or 0.1, 0.1)
    #         except:
    #             etr_today = 7.0

    #         et_today = etrf_today * etr_today
    #         et_total += et_today
    #         etr_total += etr_today

    #         is_img = "◄" if current in img_dates else ""
    #         daily_results.append({
    #             'date': date_str, 'etrf': etrf_today,
    #             'etr': etr_today, 'et': et_today})

    #         self._log(f"  {date_str}  {etrf_today:6.3f}  {etr_today:7.2f}  {et_today:7.2f}  {is_img}")

    #     self._log(f"  SUMMA: ETr={etr_total:.1f} ET={et_total:.1f} mm/oy")

    #     return daily_results, {
    #         'et_monthly': et_total,
    #         'etr_monthly': etr_total,
    #         'etrf_weighted': et_total / max(etr_total, 0.1),
    #     }
        
    # def _calc_monthly_et(self, etrf_results, geometry, start, end):
    #         """F.58: ETperiod = Σ(ETrFi × ETr24i) — PIKSEL DARAJASIDA."""
    #         start_dt = datetime.strptime(start, '%Y-%m-%d')
    #         end_dt = datetime.strptime(end, '%Y-%m-%d')
    #         n_days = (end_dt - start_dt).days + 1

    #         img_dates = [datetime.strptime(r['date'], '%Y-%m-%d') for r in etrf_results]
    #         etrf_images = [r['etrf_image'] for r in etrf_results]
    #         etrf_scalars = [r['etrf_mean'] for r in etrf_results]

    #         self._log(f"\n  ETrF nuqtalar: {len(img_dates)} ta")
    #         for d, e in zip(img_dates, etrf_scalars):
    #             self._log(f"    {d.strftime('%b %d')} → ETrF ≈ {e:.3f}")

    #         # Oylik ET xaritasi — har piksel bo'yicha
    #         et_monthly_img = ee.Image.constant(0).rename('ET_monthly')
    #         daily_results = []

    #         for day in range(n_days):
    #             current = start_dt + timedelta(days=day)
    #             date_str = current.strftime('%Y-%m-%d')

    #             # 1. ETrF IMAGE interpolatsiya (30m raster)
    #             etrf_day = self._interpolate_etrf_image(current, img_dates, etrf_images)

    #             # 2. ETr_24 IMAGE (ERA5 ~9km raster)
    #             try:
    #                 era5_day = prepare_era5_for_pm_daily(geometry, date_str)
    #                 etr_daily_img = self.ref_et.calculate(
    #                     era5_day, self.elevation, mode='daily').select('ETr')
    #             except:
    #                 etr_daily_img = ee.Image.constant(7.0).rename('ETr')

    #             # 3. ET_day = ETrF × ETr_24 (30m × 9km → 30m raster)
    #             et_day = etrf_day.multiply(etr_daily_img)

    #             # 4. Yig'indi
    #             et_monthly_img = et_monthly_img.add(et_day)

    #             # Log uchun scalar (ixtiyoriy)
    #             is_img = "◄" if current in img_dates else ""
    #             etrf_s = self._interpolate_etrf(current, img_dates, etrf_scalars)
    #             etr_s = etr_daily_img.reduceRegion(
    #                 ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
    #             ).getInfo().get('ETr', 0) or 7.0

    #             daily_results.append({
    #                 'date': date_str, 'etrf': etrf_s,
    #                 'etr': etr_s, 'et': etrf_s * etr_s})
    #             self._log(f"  {date_str}  {etrf_s:6.3f}  {etr_s:7.2f}  "
    #                     f"{etrf_s * etr_s:7.2f}  {is_img}")

    #         et_monthly_img = et_monthly_img.rename('ET_monthly')

    #         # Xarita statistikasi
    #         et_stats = et_monthly_img.reduceRegion(
    #             ee.Reducer.mean().combine(
    #                 ee.Reducer.percentile([25, 50, 75, 90]), sharedInputs=True),
    #             geometry, 300, maxPixels=1e9).getInfo()

    #         self._log(f"\n  ET xarita: mean={et_stats.get('ET_monthly_mean', 0):.1f}  "
    #                 f"P50={et_stats.get('ET_monthly_p50', 0):.1f}  "
    #                 f"P75={et_stats.get('ET_monthly_p75', 0):.1f}  "
    #                 f"P90={et_stats.get('ET_monthly_p90', 0):.1f} mm/oy")

    #         return daily_results, {
    #             'et_monthly_image': et_monthly_img,
    #             'et_monthly_mean': et_stats.get('ET_monthly_mean', 0) or 0,
    #             'et_stats': et_stats,
    #         }
    
    # def _calc_monthly_et(self, etrf_results, geometry, start, end):
    #     """F.58: ETperiod = Σ(ETrFi × ETr24i) — PIKSEL DARAJASIDA."""
    #     start_dt = datetime.strptime(start, '%Y-%m-%d')
    #     end_dt = datetime.strptime(end, '%Y-%m-%d')
    #     n_days = (end_dt - start_dt).days + 1

    #     img_dates    = [datetime.strptime(r['date'], '%Y-%m-%d') for r in etrf_results]
    #     etrf_images  = [r['etrf_image'] for r in etrf_results]
    #     etrf_scalars = [r['etrf_mean'] for r in etrf_results]
    #     crad_images  = [r['crad_image'] for r in etrf_results]
    #     crad_scalars = [r.get('crad_mean', 1.0) for r in etrf_results]

    #     # Teshik to'ldirish: composite = oylik barcha ETrF tasvirlarning piksel-wise o'rtachasi
    #     composite_etrf = (ee.ImageCollection(etrf_images)
    #                       .select('ETrF').mean().rename('ETrF'))
    #     etrf_images_filled = [img.unmask(composite_etrf) for img in etrf_images]

    #     self._log(f"\n  ETrF nuqtalar: {len(img_dates)} ta")
    #     for d, e, c in zip(img_dates, etrf_scalars, crad_scalars):
    #         self._log(f"    {d.strftime('%b %d')} → ETrF={e:.3f}  Crad={c:.3f}")

    #     # Oylik ET xaritasi — har piksel bo'yicha
    #     et_monthly_img = ee.Image.constant(0).rename('ET_monthly')
    #     daily_results = []

    #     for day in range(n_days):
    #         current = start_dt + timedelta(days=day)
    #         date_str = current.strftime('%Y-%m-%d')

    #         # 1. ETrF IMAGE interpolatsiya (30m raster) — teshiklar composite bilan to'ldirilgan
    #         etrf_day = self._interpolate_etrf_image(current, img_dates, etrf_images_filled)

    #         # 2. ETr_24 IMAGE (ERA5 ~9km raster)
    #         try:
    #             era5_day = prepare_era5_for_pm_daily(geometry, date_str)
    #             etr_daily_img = self.ref_et.calculate(
    #                 era5_day, self.elevation, mode='daily').select('ETr')
    #         except:
    #             etr_daily_img = ee.Image.constant(7.0).rename('ETr')
    #             print(f"    ⚠️ {date_str}: ETr daily image calculation failed, using fallback 7.0 mm/day")

    #         # 3. Crad — eng yaqin Landsat tasviri dan (F.56)
    #         #    Vaqt bo'yicha eng yaqin Landsat tasviri Crad ni oladi.
    #         #    Crad = f(slope, aspect, DOY) — sekin o'zgaradi, shuning uchun
    #         #    qo'shni kunlar uchun eng yaqin Landsat Crad ishlatiladi.
    #         crad_day = self._nearest_crad(current, img_dates, crad_images)

    #         # 4. F.55: ET24 = Crad × ETrF × ETr24
    #         et_day = self.et_daily.calc_et24(etrf_day, etr_daily_img, crad_day)

    #         # 4. Yig'indi
    #         et_monthly_img = et_monthly_img.add(et_day)

    #         # --- LOG TIZIMI BUTUNLAY TO'G'RILANDI ---
    #         is_img = "◄" if current in img_dates else ""
            
    #         # Xatolik beruvchi funksiyani chaqirmaymiz. Shunchaki rasterdan skalyar qiymat olamiz:
    #         try:
    #             etrf_s = etrf_day.reduceRegion(
    #                 ee.Reducer.mean(), geometry, 300, maxPixels=1e9
    #             ).getInfo().get('ETrF', 0) or 0.975
    #         except:
    #             etrf_s = 0.975
    #             print(f"    ⚠️ {date_str}: ETrF interpolation failed, using fallback 0.975")

    #         try:
    #             etr_s = etr_daily_img.reduceRegion(
    #                 ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
    #             ).getInfo().get('ETr', 0) or 7.0
    #         except:
    #             etr_s = 7.0
    #             print(f"    ⚠️ {date_str}: ETr daily image calculation failed, using fallback 7.0 mm/day")

    #         try:
    #             crad_s = crad_day.reduceRegion(
    #                 ee.Reducer.mean(), geometry, 300, maxPixels=1e9
    #             ).getInfo().get('Crad', 1.0) or 1.0
    #         except:
    #             crad_s = 1.0

    #         et_s = crad_s * etrf_s * etr_s   # F.55

    #         daily_results.append({
    #             'date': date_str,
    #             'etrf': etrf_s,
    #             'etr':  etr_s,
    #             'crad': crad_s,
    #             'et':   et_s,
    #         })
    #         self._log(f"  {date_str}  ETrF={etrf_s:.3f}  ETr={etr_s:.2f}  Crad={crad_s:.3f}  ET={et_s:.2f} mm  {is_img}")

    #     et_monthly_img = et_monthly_img.rename('ET_monthly')

    #     # Oylik o'rtacha — katta scale bilan (31 qatlamli grafik, xotira tejash)
    #     # scale=300 → GEE memory overflow; scale=5000 → 37×37 piksel, yetarli
    #     et_mean = 0.0
    #     for sc in [5000, 10000, 30000]:
    #         try:
    #             et_res = et_monthly_img.reduceRegion(
    #                 ee.Reducer.mean(), geometry, sc, maxPixels=1e9, bestEffort=True
    #             ).getInfo()
    #             v = et_res.get('ET_monthly') or et_res.get('ET_monthly_mean')
    #             if v is not None:
    #                 et_mean = float(v)
    #                 break
    #         except Exception as ex:
    #             self._log(f"  ⚠️ ET oylik mean (scale={sc}): {ex}")

    #     # Fallback: kunlik ET yig'indisidan hisoblash
    #     if et_mean == 0.0 and daily_results:
    #         et_mean = sum(r['et'] for r in daily_results)
    #         self._log(f"  ℹ️ ET oylik mean — kunlik yig'indidan: {et_mean:.1f} mm")

    #     self._log(f"\n  ET oylik mean={et_mean:.1f} mm")

    #     return daily_results, {
    #         'et_monthly_image': et_monthly_img,
    #         'et_monthly_mean':  et_mean,
    #     }
        
    def _calc_monthly_et(self, etrf_results, geometry, start, end):
        """F.58: ETperiod = Σ(ETrFi × ETr24i) — PIKSEL DARAJASIDA."""
        start_dt = datetime.strptime(start, '%Y-%m-%d')
        end_dt = datetime.strptime(end, '%Y-%m-%d')
        n_days = (end_dt - start_dt).days + 1

        img_dates    = [datetime.strptime(r['date'], '%Y-%m-%d') for r in etrf_results]
        etrf_images  = [r['etrf_image'] for r in etrf_results]
        etrf_scalars = [r['etrf_mean'] for r in etrf_results]
        crad_images  = [r['crad_image'] for r in etrf_results]
        crad_scalars = [r.get('crad_mean', 1.0) for r in etrf_results]

        # Teshik to'ldirish
        composite_etrf = (ee.ImageCollection(etrf_images)
                          .select('ETrF').mean().rename('ETrF'))
        etrf_images_filled = [img.unmask(composite_etrf) for img in etrf_images]

        # Interpolatsiya usuli
        method = getattr(self.cfg, 'etrf_interpolation', 'article')

        self._log(f"\n  ETrF nuqtalar: {len(img_dates)} ta  |  interpolatsiya: {method}")
        for d, e, c in zip(img_dates, etrf_scalars, crad_scalars):
            self._log(f"    {d.strftime('%b %d')} → ETrF={e:.3f}  Crad={c:.3f}")

        # Oylik ET xaritasi
        et_monthly_img = ee.Image.constant(0).rename('ET_monthly')
        daily_results = []

        for day in range(n_days):
            current = start_dt + timedelta(days=day)
            date_str = current.strftime('%Y-%m-%d')

            # 1. ETrF IMAGE interpolatsiya (30m raster)
            if method == 'linear':
                etrf_day = self._interpolate_etrf_image(
                    current, img_dates, etrf_images_filled)
            else:
                # 'article' — default
                etrf_day = self._interpolate_etrf_image_article(
                    current, img_dates, etrf_images_filled)

            # 2. ETr_24 IMAGE (ERA5 ~9km raster)
            try:
                era5_day = prepare_era5_for_pm_daily(geometry, date_str)
                etr_daily_img = self.ref_et.calculate(
                    era5_day, self.elevation, mode='daily').select('ETr')
            except:
                etr_daily_img = ee.Image.constant(7.0).rename('ETr')
                print(f"    ⚠️ {date_str}: ETr fallback 7.0 mm/day")

            # 3. Crad — eng yaqin Landsat tasviri dan (F.56)
            crad_day = self._nearest_crad(current, img_dates, crad_images)

            # 4. F.55: ET24 = Crad × ETrF × ETr24
            et_day = self.et_daily.calc_et24(etrf_day, etr_daily_img, crad_day)
            et_monthly_img = et_monthly_img.add(et_day)

            # --- LOG ---
            is_img = "◄" if current in img_dates else ""

            try:
                etrf_s = etrf_day.reduceRegion(
                    ee.Reducer.mean(), geometry, 300, maxPixels=1e9
                ).getInfo().get('ETrF', 0) or 0
            except:
                raise RuntimeError(f"❌ TO'XTATILDI! ETrF ni hisoblashda xatolik yuz berdi: {e}")

            try:
                etr_s = etr_daily_img.reduceRegion(
                    ee.Reducer.mean(), geometry, 5000, maxPixels=1e9
                ).getInfo().get('ETr', 0) or 7.0
            except:
                raise RuntimeError(f"❌ TO'XTATILDI! ETr ni (Referens ET) hisoblashda xatolik yuz berdi: {e}")

            try:
                crad_s = crad_day.reduceRegion(
                    ee.Reducer.mean(), geometry, 300, maxPixels=1e9
                ).getInfo().get('Crad', 1.0) or 1.0
            except:
                raise RuntimeError(f"❌ TO'XTATILDI! Crad (Bulut radiatsiyasi) ni hisoblashda xatolik yuz berdi: {e}")

            et_s = crad_s * etrf_s * etr_s

            daily_results.append({
                'date': date_str, 'etrf': etrf_s,
                'etr': etr_s, 'crad': crad_s, 'et': et_s})
            self._log(f"  {date_str}  ETrF={etrf_s:.3f}  ETr={etr_s:.2f}  "
                       f"Crad={crad_s:.3f}  ET={et_s:.2f} mm  {is_img}")

        et_monthly_img = et_monthly_img.rename('ET_monthly')

        # Oylik o'rtacha
        et_mean = 0.0
        for sc in [5000, 10000, 30000]:
            try:
                et_res = et_monthly_img.reduceRegion(
                    ee.Reducer.mean(), geometry, sc, maxPixels=1e9, bestEffort=True
                ).getInfo()
                v = et_res.get('ET_monthly') or et_res.get('ET_monthly_mean')
                if v is not None:
                    et_mean = float(v)
                    break
            except Exception as ex:
                self._log(f"  ⚠️ ET oylik mean (scale={sc}): {ex}")

        if et_mean == 0.0 and daily_results:
            et_mean = sum(r['et'] for r in daily_results)
            self._log(f"  ℹ️ ET oylik mean — kunlik yig'indidan: {et_mean:.1f} mm")

        self._log(f"\n  ET oylik mean={et_mean:.1f} mm  (interpolatsiya: {method})")

        return daily_results, {
            'et_monthly_image': et_monthly_img,
            'et_monthly_mean': et_mean,
        }


    # ================================================================
    # YORDAMCHI FUNKSIYALAR
    # ================================================================


    @staticmethod
    def _nearest_crad(target, img_dates, crad_images):
        """
        Vaqt bo'yicha eng yaqin Landsat tasviri Crad ni qaytaradi.

        Crad = f(slope, aspect, DOY) — kun sayin sekin o'zgaradi.
        Eng yaqin Landsat tasviri Crad ni ishlatish (±8 kun xatosi < 1-2%).

        Args:
            target     : datetime — hisoblash kuni
            img_dates  : list[datetime] — Landsat tasvir sanalari
            crad_images: list[ee.Image] — mos Crad tasvirlari

        Returns:
            ee.Image — 'Crad' tasvirи
        """
        if not img_dates:
            return ee.Image.constant(1.0).rename('Crad')
        if len(img_dates) == 1:
            return crad_images[0]

        min_diff = None
        nearest = crad_images[0]
        for d, img in zip(img_dates, crad_images):
            diff = abs((target - d).days)
            if min_diff is None or diff < min_diff:
                min_diff = diff
                nearest = img
        return nearest

    @staticmethod
    def _interpolate_etrf_image(target, dates, images):
        """
        ETrF IMAGE interpolatsiya — piksel bo'yicha linear.
        Birinchidan oldin/oxirgidan keyin → flat extrapolation.
        """
        if len(images) == 1:
            return images[0]

        if target <= dates[0]:
            return images[0]
        if target >= dates[-1]:
            return images[-1]

        for i in range(len(dates) - 1):
            if dates[i] <= target <= dates[i + 1]:
                total_days = (dates[i + 1] - dates[i]).days
                days_from = (target - dates[i]).days
                if total_days == 0:
                    return images[i]
                frac = days_from / total_days
                # Linear: img1 × (1-frac) + img2 × frac
                return (images[i].multiply(1 - frac)
                        .add(images[i + 1].multiply(frac)))

        return images[-1]
    
    @staticmethod
    def _interpolate_etrf_image_article(target, dates, images):
        """
        Article usuli — Allen 2007 soddalashtirilgan:
          1) Birinchi tasvirdan OLDIN → birinchi tasvir ETrF (flat)
          2) Ikki tasvir ORASIDA → (img1 + img2) / 2 (flat o'rtacha)
          3) Oxirgi tasvirdan KEYIN → oxirgi tasvir ETrF (flat)
        """
        if len(images) == 1:
            return images[0]

        # Birinchi tasvirdan oldin
        if target <= dates[0]:
            return images[0]

        # Oxirgi tasvirdan keyin
        if target >= dates[-1]:
            return images[-1]

        # Oradagi segment topish
        for i in range(len(dates) - 1):
            if dates[i] <= target <= dates[i + 1]:
                # Ikki tasvir o'rtachasi — barcha kunlar uchun BIR XIL
                return images[i].add(images[i + 1]).divide(2)

        return images[-1]

    @staticmethod
    def _interpolate_etrf_article(target, dates, etrfs):
        """Article usuli — scalar versiya (log uchun)."""
        if len(etrfs) == 1:
            return etrfs[0]
        if target <= dates[0]:
            return etrfs[0]
        if target >= dates[-1]:
            return etrfs[-1]
        for i in range(len(dates) - 1):
            if dates[i] <= target <= dates[i + 1]:
                return (etrfs[i] + etrfs[i + 1]) / 2
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

    # ================================================================
    # CLOUD PRECHECK — LULC cropland ustidagi bulut filtri
    # ================================================================

    def _cloud_precheck(self, raw_image: ee.Image, geometry) -> bool:
        """
        Rasterni ishlatish/tashlash qarorini qabul qiladi.

        Mantiq:
            1. QA_PIXEL → cloud (bit3) + cloud shadow (bit4) piksellari
            2. ESA WorldCover cropland (eq=40) mask — faqat ekin maydoni
            3. Cropland ustidagi cloud % → roi_cloud_max dan oshsa → False

        Returns:
            True  → raster yaxshi, ishlatish mumkin
            False → cropland ustida ko'p bulut, tashlasin
        """
        limit = getattr(self.cfg, 'roi_cloud_max', 30)

        try:
            # QA_PIXEL: bit3=cloud, bit4=cloud_shadow
            qa  = raw_image.select('QA_PIXEL')
            bad = (qa.bitwiseAnd(1 << 3).neq(0)
                   .Or(qa.bitwiseAnd(1 << 4).neq(0)))

            # ESA WorldCover — faqat cropland piksellar (LULC=40)
            cropland = (ee.Image("ESA/WorldCover/v100/2020")
                        .select('Map')
                        .clip(geometry)
                        .eq(40))

            # Cloud faqat cropland ustida tekshiriladi
            bad_crop = bad.updateMask(cropland)

            stats = bad_crop.reduceRegion(
                reducer    = ee.Reducer.mean(),
                geometry   = geometry,
                scale      = 120,
                maxPixels  = 1e7,
                bestEffort = True,
            ).getInfo()

            cloud_frac = stats.get('QA_PIXEL', 0) or 0
            cloud_pct  = round(cloud_frac * 100, 1)

            if cloud_pct > limit:
                self._log(f"    ⚠️  Cloud precheck FAIL: {cloud_pct}% "
                          f"(cropland, limit={limit}%) → skip")
                return False
            else:
                self._log(f"    ✅ Cloud precheck OK: {cloud_pct}% "
                          f"(cropland, limit={limit}%)")
                return True

        except Exception as e:
            self._log(f"    Cloud precheck xato ({e}) → o'tkaziladi")
            return True   # xato bo'lsa tashlama

    def _log(self, msg):
        if self.verbose:
            print(msg)
