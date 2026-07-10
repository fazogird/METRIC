"""
Quyosh geometriyasi — F.7, F.8, F.9.

F.7  cos(θrel) — qiyalik piksel uchun (slope/aspect bilan):
  cos θrel = sinδ·sinφ·coss − sinδ·cosφ·sins·cosγ
           + cosδ·cosφ·coss·cosω + cosδ·sinφ·sins·cosγ·cosω
           + cosδ·sinγ·sins·sinω

F.8  cos(θhor) — tekis yer uchun:
  cos θhor = sinδ·sinφ + cosδ·cosφ·cosω

F.9  d² = 1 / (1 + 0.033·cos(DOY·2π/365))

Aspect konversiyasi (GEE → METRIC):
  GEE:    0°=Shimol, 90°=Sharq, 180°=Janub, 270°=G'arb
  METRIC: γ=0→Janub, γ=±π→Shimol
  → cos(γ) = -cos(aspect_gee)
     sin(γ) =  sin(aspect_gee)
"""
import ee
import math


class SolarGeometry:
    """Quyosh pozitsiyasi hisoblash."""

    # def calc_solar_params(self, image_date, geometry, sun_elevation=None, dem=None):
    #     """
    #     Barcha quyosh parametrlari.

    #     Args:
    #         image_date    : str | ee.Date — tasvir sanasi
    #         geometry      : ee.Geometry   — AOI (latitude uchun)
    #         sun_elevation : ee.Number | float — Landsat metadata dan (ixtiyoriy)
    #         dem           : ee.Image | None  — 'slope_rad', 'aspect_rad' bandlari
    #                         None → F.8 (tekis yer); berilsa → F.7 (qiyalik)

    #     Returns:
    #         dict:
    #           doy           — ee.Number
    #           d2            — ee.Number  (F.9)
    #           delta         — ee.Number  (rad, deklinatsiya)
    #           omega         — ee.Number  (rad, soat burchagi ~10:30)
    #           cos_theta_hor — ee.Number  (F.8, tekis, skalyar)
    #           cos_theta_rel — ee.Number | ee.Image
    #                           dem=None → F.8 skalyar (tekis yer)
    #                           dem≠None → F.7 ee.Image (har piksel uchun)
    #           sun_elevation — ee.Number  (daraja)
    #     """
    #     date = ee.Date(image_date)
    #     doy  = date.getRelative('day', 'year').add(1)

    #     # ── F.9: d² ───────────────────────────────────────────────
    #     d2 = (doy.multiply(2 * math.pi / 365).cos()
    #           .multiply(0.033).add(1).pow(-1))

    #     # ── Deklinatsiya (Cooper tenglamasi) ──────────────────────
    #     delta = (doy.subtract(1).multiply(2 * math.pi / 365)
    #              .subtract(ee.Number(math.pi * 80 / 365))
    #              .sin().multiply(0.4093))

    #     # ── Soat burchagi — Landsat overpass ~10:30 mahalliy ─────
    #     # ω ≈ -0.3 rad  (tushga 22 daqiqa qolgan)
    #     omega = ee.Number(-0.3)

    #     # ── F.8: cos(θhor) — tekis yer (skalyar) ─────────────────
    #     lat_deg = ee.Number(geometry.centroid().coordinates().get(1))
    #     lat_rad = lat_deg.multiply(math.pi / 180)

    #     cos_theta_hor = (delta.sin().multiply(lat_rad.sin())
    #                      .add(delta.cos().multiply(lat_rad.cos())
    #                           .multiply(omega.cos())))

    #     # ── Sun elevation ──────────────────────────────────────────
    #     if sun_elevation is not None:
    #         sun_elev = ee.Number(sun_elevation)
    #     else:
    #         # cos(θ) dan sun elevation hisoblash
    #         sun_elev = (cos_theta_hor.acos()
    #                     .multiply(-1).add(math.pi / 2)
    #                     .multiply(180 / math.pi))

    #     # ── cos(θrel): F.7 yoki F.8 ───────────────────────────────
    #     if dem is not None:
    #         cos_theta_rel = self._calc_cos_theta_rel_f7(dem, delta, omega)
    #     else:
    #         # DEM yo'q → tekis yer (F.8), skalyar
    #         cos_theta_rel = cos_theta_hor

    #     return {
    #         'doy':           doy,
    #         'd2':            d2,
    #         'delta':         delta,
    #         'omega':         omega,
    #         'cos_theta_hor': cos_theta_hor,   # ee.Number — tekis (F.8)
    #         'cos_theta_rel': cos_theta_rel,   # ee.Image (F.7) yoki ee.Number (F.8)
    #         'sun_elevation': sun_elev,
    #     }
        
    def calc_solar_params(self, date_str, geometry, sun_elevation=None,
                          dem=None, img_date=None, timezone_lon_deg=0.0):
        """
        Quyosh geometriyasi parametrlari.

        Args:
            date_str:  'YYYY-MM-dd'
            geometry:  ee.Geometry
            sun_elevation: ee.Number (Landsat metadata dan)
            dem:       ee.Image (slope_rad, aspect_rad bandlari)
            img_date:  ee.Date — Landsat to'liq vaqti (soat bilan)
            timezone_lon_deg: float — vaqt zonasi markaz meridian
                              Idaho: -105, O'zbekiston: +75
        """
        # ── DOY ───────────────────────────────────────────────────
        doy = ee.Date(date_str).getRelative('day', 'year').add(1)

        # ── d² (yer-quyosh masofa koeffitsienti) ─────────────────
        d2 = (doy.multiply(2 * math.pi / 365)
              .cos().multiply(0.033).add(1).pow(-1))

        # ── Deklinatsiya ──────────────────────────────────────────
        delta = (doy.subtract(1).multiply(2 * math.pi / 365)
                 .subtract(ee.Number(math.pi * 80 / 365))
                 .sin().multiply(0.4093))

        # ── Omega — piksel bo'yicha (RefETCalculator usulida) ─────
        if img_date is not None:
            # Landsat vaqtidan soat olish (UTC)
            hour_utc = ee.Number(img_date.get('hour'))

            # Mavsumiy tuzatma (Sc)
            b_sc = doy.subtract(81).multiply(2 * math.pi / 364)
            Sc = (b_sc.multiply(2).sin().multiply(0.1645)
                  .subtract(b_sc.cos().multiply(0.1255))
                  .subtract(b_sc.sin().multiply(0.025)))

            # Piksel uzunligi (ee.Image)
            lon_img = ee.Image.pixelLonLat().select('longitude')
            Lz = ee.Number(timezone_lon_deg)

            # Solar time = hour_UTC + (Lm - Lz)/15 + Sc
            # (0.06667 = 1/15)
            sol_time = (lon_img.subtract(Lz)
                        .multiply(0.06667)
                        .add(hour_utc)
                        .add(Sc))

            # omega = (sol_time - 12) × π/12
            omega = sol_time.subtract(12).multiply(math.pi / 12)
            # omega endi ee.Image — har piksel o'z qiymati!
        else:
            # Fallback — Landsat 10:30 local
            omega = ee.Number(-0.393)
            print(f"    ⚠️ {date_str}: img_date berilmagan, omega fallback -0.393 rad ishlatiladi")

        # ── cos(θhor) ─────────────────────────────────────────────
        lat_rad = (ee.Image.pixelLonLat().select('latitude')
                   .multiply(math.pi / 180))

        # omega ee.Image bo'lsa, cos_theta_hor ham ee.Image bo'ladi
        cos_theta_hor = (lat_rad.sin().multiply(delta.sin())
                         .add(lat_rad.cos().multiply(delta.cos())
                              .multiply(omega.cos())))

        # ── Sun elevation ─────────────────────────────────────────
        if sun_elevation is not None:
            sun_elev = ee.Number(sun_elevation)
        else:
            # cos(θ) dan hisoblash — o'rtacha qiymat olish kerak
            sun_elev_img = (cos_theta_hor.acos()
                            .multiply(-1).add(math.pi / 2)
                            .multiply(180 / math.pi))
            sun_elev = ee.Number(sun_elev_img.reduceRegion(
                ee.Reducer.mean(), geometry, 1000, maxPixels=1e9
            ).values().get(0))

        # ── cos(θrel): F.7 yoki F.8 ──────────────────────────────
        if dem is not None:
            cos_theta_rel = self._calc_cos_theta_rel_f7(dem, delta, omega)
        else:
            cos_theta_rel = cos_theta_hor

        return {
            'doy':           doy,
            'd2':            d2,
            'delta':         delta,
            'omega':         omega,           # ee.Image (piksel-wise)
            'cos_theta_hor': cos_theta_hor,   # ee.Image (piksel-wise)
            'cos_theta_rel': cos_theta_rel,   # ee.Image
            'sun_elevation': sun_elev,        # ee.Number
        }


    # ================================================================
    # F.7 — Qiyalik piksel uchun cos(θrel)
    # ================================================================

    @staticmethod
    def _calc_cos_theta_rel_f7(dem, delta, omega):
        """
        F.7: cos(θrel) — har piksel uchun slope/aspect tuzatmasi.

        cos θrel = sinδ·sinφ·coss − sinδ·cosφ·sins·cosγ
                 + cosδ·cosφ·coss·cosω + cosδ·sinφ·sins·cosγ·cosω
                 + cosδ·sinγ·sins·sinω

        Args:
            dem   : ee.Image — 'slope_rad', 'aspect_rad' bandlari
            delta : ee.Number — deklinatsiya (rad)
            omega : ee.Number — soat burchagi (rad)

        Returns:
            ee.Image — cos(θrel), manfiy = soya
        """
        lat_rad = ee.Image.pixelLonLat().select('latitude').multiply(math.pi / 180)

        slope_s = dem.select('slope_rad')      # s (rad)
        asp_rad = dem.select('aspect_rad')      # GEE: 0=N, 90=E, 180=S, 270=W (rad)

        # GEE → METRIC γ:  cos(γ) = -cos(asp),  sin(γ) = sin(asp)
        cos_gamma = asp_rad.cos().multiply(-1)
        sin_gamma = asp_rad.sin()

        sin_d = delta.sin()
        cos_d = delta.cos()
        sin_p = lat_rad.sin()
        cos_p = lat_rad.cos()
        cos_s = slope_s.cos()
        sin_s = slope_s.sin()
        cos_w = omega.cos()
        sin_w = omega.sin()

        # ⚠️ GEE Python: Number.multiply(Image) → Number (noto'g'ri!)
        #               Image.multiply(Number) → Image (to'g'ri!)
        # → har bir hadda Image (sin_p/cos_p/cos_s/sin_s) CHAPDA bo'lishi kerak
        t1 =  sin_p.multiply(sin_d).multiply(cos_s)
        t2 =  cos_p.multiply(sin_d).multiply(sin_s).multiply(cos_gamma).multiply(-1)
        t3 =  cos_p.multiply(cos_d).multiply(cos_s).multiply(cos_w)
        t4 =  sin_p.multiply(cos_d).multiply(sin_s).multiply(cos_gamma).multiply(cos_w)
        t5 =  sin_s.multiply(cos_d).multiply(sin_gamma).multiply(sin_w)

        return t1.add(t2).add(t3).add(t4).add(t5).rename('cos_theta_rel')
