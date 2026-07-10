"""
core/reference_et.py — ASCE Standardized Reference ET.

ASCE-EWRI 2005 to'liq implementatsiya.
calculate() har doim ikkala band qaytaradi:
  ETr — alfalfa (tall) reference
  ETo — grass   (short) reference

Kirish bandlari: T_air(°C), P_kPa, u2(m/s), ea(kPa), Rs(MJ/m²/hr yoki /day)
"""
import ee
import math


class RefETCalculator:
    """ASCE Standardized Penman-Monteith — ETr (alfalfa) va ETo (grass)."""

    Gsc    = 0.0820   # MJ m⁻² min⁻¹
    albedo = 0.23

    # ASCE coefficients — alfalfa (tall)  [ASCE-EWRI 2005, Table 1]
    _ETR = dict(Cn_hr=66.0,  Cd_hr_day=0.25, Cd_hr_night=1.7,
                G_hr_day=0.04, G_hr_night=0.20,
                Cn_day=1600.0, Cd_day=0.38)

    # ASCE coefficients — grass (short)
    _ETO = dict(Cn_hr=37.0,  Cd_hr_day=0.24, Cd_hr_night=0.96,
                G_hr_day=0.10, G_hr_night=0.50,
                Cn_day=900.0,  Cd_day=0.34)

    def __init__(self, timezone_lon_deg=0.0, ref_type='alfalfa'):
        self.Lz       = ee.Number(timezone_lon_deg)
        self.ref_type = ref_type   # tashqaridan so'ralganda ishlatiladi

    # ================================================================
    # YORDAMCHI METODLAR
    # ================================================================

    @staticmethod
    def _esat(T_c):
        """Saturation vapor pressure [kPa] — ASCE Eq. 7."""
        return T_c.expression('0.6108 * exp((17.27 * T) / (T + 237.3))', {'T': T_c})

    @staticmethod
    def _delta_svp(T_c, es_kPa):
        """Slope of SVP curve [kPa/°C] — ASCE Eq. 5."""
        return es_kPa.multiply(4098).divide(T_c.add(237.3).pow(2))

    @staticmethod
    def _gamma(P_kPa):
        """Psychrometric constant [kPa/°C] — ASCE Eq. 4."""
        return P_kPa.multiply(0.000665)

    @staticmethod
    def _doy_hour(img):
        img  = ee.Image(img)
        date = ee.Date(img.get('system:time_start'))
        doy  = ee.Number(date.getRelative('day', 'year')).add(1)
        hour = ee.Number(date.get('hour'))
        return doy, hour

    @staticmethod
    def _dr_decl(doy):
        """dr (inverse earth-sun distance) va delta (declination)."""
        J   = ee.Number(doy)
        dr  = J.multiply(2 * math.pi / 365.0).cos().multiply(0.033).add(1.0)
        dec = J.multiply(2 * math.pi / 365.0).subtract(1.39).sin().multiply(0.409)
        return dr, dec

    @staticmethod
    def _seasonal_correction(doy):
        """Seasonal correction Sc [hr] — ASCE Eq. 58."""
        J  = ee.Number(doy)
        b  = J.subtract(81).multiply(2 * math.pi / 364.0)
        Sc = (b.multiply(2).sin().multiply(0.1645)
              .subtract(b.cos().multiply(0.1255))
              .subtract(b.sin().multiply(0.025)))
        return Sc

    def _omega(self, lon_deg, doy, hour, offset=0.5):
        """Soat burchagi omega [rad] — ASCE Eq. 55.

        offset=0.5  → davr o'rtasi (Ra uchun; RefET hourly.py:91)
        offset=0.0  → davr boshi (fcd beta uchun; RefET calcs.py:155-156)
        """
        Sc       = self._seasonal_correction(doy)
        t_local  = ee.Number(hour).add(self.Lz.divide(15)).add(offset)
        lon_img  = ee.Image(lon_deg)
        sol_time = (lon_img.subtract(ee.Image.constant(self.Lz))
                    .multiply(0.06667)
                    .add(t_local)
                    .add(Sc))
        omega = sol_time.subtract(12).multiply(math.pi / 12.0)
        # RefET calcs.py:348 — wrap to [-π, π]
        return omega.add(math.pi).mod(2 * math.pi).subtract(math.pi)

    def _omega_mid_hour(self, lon_deg, doy, hour):
        """Davr o'rtasidagi omega — Ra hisoblash uchun."""
        return self._omega(lon_deg, doy, hour, offset=0.5)

    # ================================================================
    # RADIATSIYA
    # ================================================================

    def Ra_hourly(self, lat_rad, lon_deg, doy, hour, dr, dec):
        """Hourly extraterrestrial radiation [MJ m⁻² h⁻¹] — ASCE Eq. 48."""
        # _omega_mid_hour ichida wrap allaqachon qo'llanilgan
        omega  = self._omega_mid_hour(lon_deg, doy, hour)
        omega1 = omega.subtract(math.pi / 24.0)
        omega2 = omega.add(math.pi / 24.0)
        # RefET calcs.py:387 — clip argument to [-1,1] before arccos (polar safety)
        ws     = lat_rad.tan().multiply(-1).multiply(dec.tan()).clamp(-1, 1).acos()
        omega1c = omega1.max(ws.multiply(-1)).min(ws)
        omega2c = omega2.max(ws.multiply(-1)).min(ws)
        term1  = omega2c.subtract(omega1c).multiply(lat_rad.sin()).multiply(dec.sin())
        term2  = lat_rad.cos().multiply(dec.cos()).multiply(omega2c.sin().subtract(omega1c.sin()))
        Ra     = term1.add(term2).multiply(dr).multiply(12 * 60 / math.pi * self.Gsc)
        return Ra.max(0)

    def Ra_daily(self, lat_rad, doy, dr, dec):
        """Daily extraterrestrial radiation [MJ m⁻² d⁻¹] — ASCE Eq. 21."""
        # RefET calcs.py:387 — clip argument to [-1,1] before arccos (polar safety)
        ws    = lat_rad.tan().multiply(-1).multiply(dec.tan()).clamp(-1, 1).acos()
        term1 = ws.multiply(lat_rad.sin()).multiply(dec.sin())
        term2 = lat_rad.cos().multiply(dec.cos()).multiply(ws.sin())
        Ra    = term1.add(term2).multiply(dr).multiply(24 * 60 / math.pi * self.Gsc)
        return Ra

    @staticmethod
    def Rso(Ra, z_m):
        """Clear-sky solar radiation — ASCE Eq. 19/45: (0.75 + 2e-5*z)*Ra."""
        return z_m.multiply(2e-5).add(0.75).multiply(Ra)

    @staticmethod
    def Rnl_hourly(T_k, ea_kPa, Rs, Rso, lat_rad=None, dec=None, omega=None):
        """Hourly net long-wave radiation — ASCE Eq. 44.

        T_k      : T [°C] + 273.16  (RefET calcs.py:730)
        lat_rad  : ee.Image — latitude [rad] (beta hisoblash uchun)
        dec      : ee.Number — deklinatsiya [rad]
        omega    : ee.Image — soat burchagi DAVR BOSHIDA [rad] (RefET fcd uchun)
        """
        Rso_safe = Rso.max(1e-6)
        rs_rso   = Rs.divide(Rso_safe).clamp(0.3, 1.0)
        fcd      = rs_rso.multiply(1.35).subtract(0.35)

        # RefET calcs.py:669-681 — beta < 0.3 rad da fcd=1 (quyosh past)
        if lat_rad is not None and dec is not None and omega is not None:
            sin_beta = (lat_rad.sin().multiply(dec.sin())
                        .add(lat_rad.cos().multiply(dec.cos()).multiply(omega.cos())))
            beta = sin_beta.clamp(-1, 1).asin()
            fcd  = fcd.where(beta.lt(0.3), ee.Image.constant(1.0))
        else:
            # Fallback: Rs ≤ 0.01 → tun (lat/omega berilmagan holda)
            fcd  = fcd.where(Rs.lte(0.01), ee.Image.constant(1.0))

        emiss = ea_kPa.sqrt().multiply(-0.14).add(0.34)
        # RefET calcs.py:730 — T + 273.16 (T_k = T + 273.16 keladi)
        return T_k.pow(4).multiply(2.042e-10).multiply(emiss).multiply(fcd)

    @staticmethod
    def Rnl_daily(Tmax_k, Tmin_k, ea_kPa, Rs, Rso):
        """Daily net long-wave radiation — ASCE Eq. 17."""
        # RefET: kunlik ham [0.3, 1.0] — soatlik bilan bir xil chegaralar
        rs_rso = Rs.divide(Rso.max(1e-6)).clamp(0.3, 1.0)
        fcd    = rs_rso.multiply(1.35).subtract(0.35)
        emiss  = ea_kPa.sqrt().multiply(-0.14).add(0.34)
        tmean4 = Tmax_k.pow(4).add(Tmin_k.pow(4)).divide(2.0)
        # σ = 4.901e-9 MJ m⁻² day⁻¹ K⁻⁴  (RefET: 4.901, ASCE rounded 4.903)
        return ee.Image.constant(4.901e-9).multiply(tmean4).multiply(emiss).multiply(fcd)

    @staticmethod
    def _soil_heat_flux(Rn, mode, coef):
        """G = Rn × koeffitsient (kun/tun). Kunlik uchun G=0."""
        if mode == 'daily':
            return ee.Image.constant(0)
        g_day   = Rn.multiply(coef['G_hr_day'])
        g_night = Rn.multiply(coef['G_hr_night'])
        return g_day.where(Rn.lt(0), g_night)

    @staticmethod
    def _pm_hourly(Delta, Gamma, Rn, G, u2, vpd, T, coef, is_night):
        """Penman-Monteith soatlik — kun/tun alohida Cd bilan."""
        Cn      = coef['Cn_hr']
        num1    = Delta.multiply(Rn.subtract(G)).multiply(0.408)
        num2    = Gamma.multiply(u2).multiply(vpd).multiply(
                      ee.Image.constant(Cn).divide(T.add(273)))
        et_day  = num1.add(num2).divide(
                      Delta.add(Gamma.multiply(u2.multiply(coef['Cd_hr_day']).add(1.0))))
        et_night = num1.add(num2).divide(
                      Delta.add(Gamma.multiply(u2.multiply(coef['Cd_hr_night']).add(1.0))))
        return et_day.where(is_night, et_night).max(0)

    @staticmethod
    def _pm_daily(Delta, Gamma, Rn, u2, vpd, T, coef):
        """Penman-Monteith kunlik — ASCE Eq. 1."""
        num1 = Delta.multiply(Rn).multiply(0.408)
        num2 = Gamma.multiply(u2).multiply(vpd).multiply(
                   ee.Image.constant(coef['Cn_day']).divide(T.add(273)))
        den  = Delta.add(Gamma.multiply(u2.multiply(coef['Cd_day']).add(1.0)))
        return num1.add(num2).divide(den).max(0)

    # ================================================================
    # ASOSIY HISOBLASH
    # ================================================================

    def calculate(self, img, dem, mode='hourly'):
        """
        ETr (alfalfa) va ETo (grass) — bir o'tishda hisoblash.

        Args:
            img  : ee.Image — T_air(°C), P_kPa, u2(m/s), ea(kPa), Rs(MJ/m²)
                              kunlik bo'lsa + T_max, T_min (ixtiyoriy)
            dem  : ee.Image — elevation [m]
            mode : 'hourly' | 'daily'

        Returns:
            ee.Image — kirish bandlari + ETr + ETo
                       (+ Ra, Rso, Rnl, Rn_ref, G_ref oraliq bandlari)
        """
        img = ee.Image(img)
        dem = ee.Image(dem)

        ll      = ee.Image.pixelLonLat()
        lat_rad = ll.select('latitude').multiply(math.pi / 180.0)
        lon_deg = ll.select('longitude')

        doy, hour = self._doy_hour(img)
        dr,  dec  = self._dr_decl(doy)

        # Extraterrestrial radiation
        if mode == 'hourly':
            Ra = self.Ra_hourly(lat_rad, lon_deg, doy, hour, dr, dec)
        else:
            Ra = self.Ra_daily(lat_rad, doy, dr, dec)

        Rso_img = self.Rso(Ra, dem)

        T  = img.select('T_air')
        u2 = img.select('u2')
        P  = img.select('P_kPa')
        ea = img.select('ea')
        Rs = img.select('Rs')

        # Saturation vapor pressure
        if mode == 'daily':
            bnames  = img.bandNames()
            T_max_c = ee.Image(ee.Algorithms.If(bnames.contains('T_max'), img.select('T_max'), T))
            T_min_c = ee.Image(ee.Algorithms.If(bnames.contains('T_min'), img.select('T_min'), T))
            es = self._esat(T_max_c).add(self._esat(T_min_c)).divide(2).rename('es')
        else:
            es = self._esat(T).rename('es')

        # RefET: Delta Tmean dan to'g'ridan hisoblanadi (avg(esat) emas — Eq.5)
        es_for_delta = self._esat(T)
        Delta = self._delta_svp(T, es_for_delta).rename('Delta')
        Gamma = self._gamma(P).rename('Gamma')

        # RefET: soatlik — clamp YO'Q (hourly.py:135); kunlik — max(0) (daily.py:152)
        if mode == 'hourly':
            vpd = es.subtract(ea)
        else:
            vpd = es.subtract(ea).max(0)

        # Net radiation (alfalfa va grass uchun bir xil: albedo=0.23)
        Rns = Rs.multiply(1.0 - self.albedo)
        if mode == 'hourly':
            # omega davr BOSHI (fcd beta uchun; RefET calcs.py:155-156)
            omega_start = self._omega(lon_deg, doy, hour, offset=0.0)
            # RefET calcs.py:730 — T + 273.16
            Rnl = self.Rnl_hourly(
                T.add(273.16), ea, Rs, Rso_img,
                lat_rad=lat_rad, dec=dec, omega=omega_start)
        else:
            bnames = img.bandNames()
            # RefET calcs.py:708 — T + 273.16
            Tmax_k = ee.Image(ee.Algorithms.If(
                bnames.contains('T_max'), img.select('T_max').add(273.16), T.add(273.16)))
            Tmin_k = ee.Image(ee.Algorithms.If(
                bnames.contains('T_min'), img.select('T_min').add(273.16), T.add(273.16)))
            Rnl = self.Rnl_daily(Tmax_k, Tmin_k, ea, Rs, Rso_img)

        Rn       = Rns.subtract(Rnl).rename('Rn_ref')
        is_night = Rn.lt(0)

        # ── ETr — alfalfa (tall) ──────────────────────────────────
        G_etr = self._soil_heat_flux(Rn, mode, self._ETR)
        if mode == 'hourly':
            ETr = self._pm_hourly(Delta, Gamma, Rn, G_etr, u2, vpd, T, self._ETR, is_night)
        else:
            ETr = self._pm_daily(Delta, Gamma, Rn, u2, vpd, T, self._ETR)
        ETr = ETr.rename('ETr')

        # ── ETo — grass (short) ───────────────────────────────────
        G_eto = self._soil_heat_flux(Rn, mode, self._ETO)
        if mode == 'hourly':
            ETo = self._pm_hourly(Delta, Gamma, Rn, G_eto, u2, vpd, T, self._ETO, is_night)
        else:
            ETo = self._pm_daily(Delta, Gamma, Rn, u2, vpd, T, self._ETO)
        ETo = ETo.rename('ETo')

        return img.addBands([es, Delta, Gamma,
                             Ra.rename('Ra'), Rso_img.rename('Rso'),
                             Rnl.rename('Rnl'), Rn,
                             G_etr.rename('G_ref'),
                             ETr, ETo])
